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
    CUvideodecoder cudecoder;
    CUvideoparser cuparser;
    int num_surfaces;
} nvdec_ctx_t;

static int CUDAAPI HandleVideoSequence(void *opaque, CUVIDEOFORMAT *format)
{
    decoder_t *p_dec = (decoder_t *) opaque;
    nvdec_ctx_t *p_ctx = p_dec->p_sys;

    // post processed output's dimensions need to be 2-aligned in NVDEC
    p_dec->fmt_out.video.i_width = (format->coded_width + 1) & ~1;
    p_dec->fmt_out.video.i_height = (format->coded_height + 1) & ~1;

    CUVIDDECODECREATEINFO dparams = {0};
    dparams.ulWidth             = format->coded_width;
    dparams.ulHeight            = format->coded_height;
    dparams.ulTargetWidth       = p_dec->fmt_out.video.i_width;
    dparams.ulTargetHeight      = p_dec->fmt_out.video.i_height;
    dparams.bitDepthMinus8      = format->bit_depth_luma_minus8;
    dparams.OutputFormat        = cudaVideoSurfaceFormat_NV12;
    dparams.CodecType           = format->codec;
    dparams.ChromaFormat        = format->chroma_format;
    dparams.ulNumDecodeSurfaces = p_ctx->num_surfaces;
    dparams.ulNumOutputSurfaces = 1;

    CUresult result = cuvidCreateDecoder(&p_ctx->cudecoder, &dparams);

    if (result != CUDA_SUCCESS) {
        msg_Err(p_dec, "Could not create decoder object");
        return 0;
    }
    return 1;
}

static int CUDAAPI HandlePictureDecode(void *opaque, CUVIDPICPARAMS *params)
{
    decoder_t *p_dec = (decoder_t *) opaque;
    nvdec_ctx_t *p_ctx = p_dec->p_sys;

    CUresult result = cuvidDecodePicture(&p_ctx->cudecoder, params);

    if (result != CUDA_SUCCESS) {
        msg_Err(p_dec, "Could not create decoder object");
        return 0;
    }
    return 1;
}

static int CUDAAPI HandlePictureDisplay(void *opaque, CUVIDPARSERDISPINFO *dispinfo)
{
    return 1;
}

static int DecodeBlock(decoder_t *p_dec, block_t *p_block)
{
    nvdec_ctx_t *p_ctx = p_dec->p_sys;
    if (p_block == NULL) {
        // TODO flush
        return VLCDEC_SUCCESS;
    }

    CUVIDSOURCEDATAPACKET cupacket = {0};
    cupacket.payload_size = p_block->i_buffer;
    cupacket.payload = p_block->p_buffer;
    cupacket.timestamp = p_block->i_pts;

    CUresult result = cuvidParseVideoData(&p_ctx->cuparser, &cupacket);

    if (result != CUDA_SUCCESS) {
        msg_Err(p_dec, "Could not send packet to NVDEC parser");
        return VLCDEC_ECRITICAL;
    } else {
        msg_Err(p_dec, "Sent packet to NVDEC parser");
    }
    return VLCDEC_SUCCESS;
}

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
    if (result != CUDA_SUCCESS || !caps.bIsSupported) {
        msg_Err(p_dec, "No hardware for NVDEC");
        return VLC_EGENERIC;
    }

    p_ctx->num_surfaces = MAX_SURFACES;

    CUVIDPARSERPARAMS pparams = {0};
    pparams.CodecType               = cudaVideoCodec_H264;
    pparams.ulClockRate             = 1e6;
    pparams.ulMaxDisplayDelay       = 4;
    pparams.ulMaxNumDecodeSurfaces  = p_ctx->num_surfaces;
    pparams.pUserData               = p_dec;
    pparams.pfnSequenceCallback     = HandleVideoSequence;
    pparams.pfnDecodePicture        = HandlePictureDecode;
    pparams.pfnDisplayPicture       = HandlePictureDisplay;
    result = cuvidCreateVideoParser(&p_ctx->cuparser, &pparams);
    if (result != CUDA_SUCCESS) {
        msg_Err(p_dec, "Could not create parser object");
        return VLC_EGENERIC;
    }
    return VLC_SUCCESS;
}

static void CloseDecoder(vlc_object_t *p_this)
{
    decoder_t *p_dec = (decoder_t *) p_this;
    nvdec_ctx_t *p_ctx = p_dec->p_sys;
    cuvidDestroyDecoder(&p_ctx->cudecoder);
    cuvidDestroyVideoParser(&p_ctx->cuparser);
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
