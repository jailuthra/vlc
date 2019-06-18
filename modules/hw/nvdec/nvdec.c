/*****************************************************************************
 * nvdec.c: NVDEC hw video decoder
 *****************************************************************************
 * Copyright (C) 2019 VLC authors and VideoLAN
 *
 * Authors: Jai Luthra <me@jailuthra.in>
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston MA 02110-1301, USA.
 *****************************************************************************/
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <vlc_plugin.h>
#include <vlc_common.h>
#include <vlc_codec.h>
#include <vlc_messages.h>
#include "nvcuvid.h"

#define MAX_SURFACES 25

typedef struct nvdec_ctx {
    CUvideodecoder decoder;
    CUstream stream;
    int num_surfaces;
} nvdec_ctx_t;

static int DecodeBlock(decoder_t *p_dec, block_t *p_block);

static int OpenDecoder(vlc_object_t *p_this)
{
    decoder_t *p_dec = (decoder_t *) p_this;
    nvdec_ctx_t *p_ctx;
    p_ctx = calloc(1, sizeof(*p_ctx));
    if (!p_ctx) {
        return VLC_ENOMEM;
    }
    p_dec->p_sys = p_ctx;
    p_dec->pf_decode = DecodeBlock;

    if (p_dec->fmt_in.i_codec != VLC_CODEC_H264)
        return VLC_EGENERIC;

    CUVIDDECODECAPS caps = {0};
    caps.eCodecType         = cudaVideoCodec_H264;
    caps.eChromaFormat      = cudaVideoSurfaceFormat_NV12;
    caps.nBitDepthMinus8    = p_dec->fmt_in.video.i_bits_per_pixel - 8;
    CUresult result = cuvidGetDecoderCaps(&caps);
    if (result != CUDA_SUCCESS && caps.bIsSupported) {
        msg_Err(p_dec, "No hardware for NVDEC");
        return VLC_EGENERIC;
    } else {
        msg_Dbg(p_dec, "Hardware is capable of decoding");
    }
    
    p_ctx->num_surfaces = MAX_SURFACES;
    CUVIDDECODECREATEINFO params = {0};
    params.ulWidth             = p_dec->fmt_in.video.i_width;
    params.ulHeight            = p_dec->fmt_in.video.i_height;
    params.ulTargetWidth       = p_dec->fmt_in.video.i_width;
    params.ulTargetHeight      = p_dec->fmt_in.video.i_height;
    params.bitDepthMinus8      = p_dec->fmt_in.video.i_bits_per_pixel - 8;
    params.OutputFormat        = cudaVideoSurfaceFormat_NV12;
    params.CodecType           = cudaVideoCodec_H264;
    params.ChromaFormat        = cudaVideoChromaFormat_420;
    params.ulNumDecodeSurfaces = p_ctx->num_surfaces;
    params.ulNumOutputSurfaces = 1;

    result = cuvidCreateDecoder(&p_ctx->decoder, &params);
    if (result != CUDA_SUCCESS) {
        msg_Err(p_dec, "Could not create decoder object");
        return VLC_EGENERIC;
    } else {
        msg_Dbg(p_dec, "Created decoder object");
    }

    return VLC_SUCCESS;
}

static void CloseDecoder(vlc_object_t *p_this)
{
    decoder_t *p_dec = (decoder_t *) p_this;
    nvdec_ctx_t *p_ctx = p_dec->p_sys;
    cuvidDestroyDecoder(&p_ctx->decoder);
    free(p_ctx);
}

vlc_module_begin ()
    set_description(N_("NVDEC video decoder"))
    set_shortname("nvdec")
    set_capability("video decoder", 700)
    set_category(CAT_INPUT)
    set_subcategory(SUBCAT_INPUT_VCODEC)
    set_callbacks(OpenDecoder, CloseDecoder)
vlc_module_end ()
