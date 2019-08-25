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

#include <vlc_common.h>
#include <vlc_plugin.h>
#include <vlc_codec.h>
#include <vlc_messages.h>
#include <ffnvcodec/dynlink_cuda.h>
#include <ffnvcodec/dynlink_loader.h>
#include "hxxx_helper.h"

#define NVDEC_PIPELINE_SIZE 4
#define MAX_HXXX_SURFACES (18 + 2)

typedef struct nvdec_ctx {
    CuvidFunctions              *functions;
    CudaFunctions               *cudaFunctions;
    CUcontext                   cuCtx;
    CUvideodecoder              cudecoder;
    CUvideoparser               cuparser;
    CUVIDPARSERPARAMS           pparams;
    struct hxxx_helper          hh;
    bool                        b_is_hxxx;
    int                         i_nb_surface; ///<  number of GPU surfaces allocated
    bool                        b_xps_pushed; ///< (for xvcC) parameter sets pushed (SPS/PPS/VPS)
    bool                        b_progressive;
    cudaVideoDeinterlaceMode    deintMode;
} nvdec_ctx_t;

static int OpenDecoder(vlc_object_t *p_this);
static void CloseDecoder(vlc_object_t *p_this);

vlc_module_begin ()
    set_description(N_("NVDEC video decoder"))
    set_shortname("nvdec")
    set_capability("video decoder", 1000)
    set_category(CAT_INPUT)
    set_subcategory(SUBCAT_INPUT_VCODEC)
    set_callbacks(OpenDecoder, CloseDecoder)
vlc_module_end ()

static inline int CudaCall(decoder_t *p_dec, CUresult result, const char *psz_func)
{
    if (unlikely(result != CUDA_SUCCESS)) {
        const char *psz_err;
        nvdec_ctx_t *p_ctx = p_dec->p_sys;
        p_ctx->cudaFunctions->cuGetErrorName(result, &psz_err);
        msg_Err(p_dec, "%s: %s failed", psz_err, psz_func);
        return VLC_EGENERIC;
    }
    return VLC_SUCCESS;
}

#define CUDA_CALL(func) CudaCall(p_dec, (func), #func)

static int UpdateOutputFormat(decoder_t *p_dec, CUVIDEOFORMAT *p_format)
{
    // post processed output's dimensions need to be 2-aligned in NVDEC
    p_dec->fmt_out.video.i_width = (p_format->coded_width + 1) & ~1;
    p_dec->fmt_out.video.i_height = (p_format->coded_height + 1) & ~1;
    // frame rate
    p_dec->fmt_out.video.i_frame_rate = p_format->frame_rate.numerator;
    p_dec->fmt_out.video.i_frame_rate_base = p_format->frame_rate.denominator;
    // bit depth and chroma
    unsigned int i_bpp = p_format->bit_depth_luma_minus8 + 8;
    vlc_fourcc_t i_chroma;
    if (i_bpp == 8) {
        switch (p_format->chroma_format) {
            case cudaVideoChromaFormat_420:
                i_chroma = VLC_CODEC_NV12;
                break;
            case cudaVideoChromaFormat_444:
                i_chroma = VLC_CODEC_I444;
                break;
            default:
                return VLC_EGENERIC;
        }
    } else {
        // TODO: figure out how to support P016 or YUV444_16bit
        return VLC_EGENERIC;
    }
    p_dec->fmt_out.video.i_bits_per_pixel = i_bpp;
    p_dec->fmt_out.video.i_chroma = p_dec->fmt_out.i_codec = i_chroma;
    decoder_UpdateVideoFormat(p_dec);
    return VLC_SUCCESS;
}

static int MapSurfaceFmt(int i_vlc_fourcc)
{
    switch (i_vlc_fourcc) {
        case VLC_CODEC_NV12: return cudaVideoSurfaceFormat_NV12;
        case VLC_CODEC_I444: return cudaVideoSurfaceFormat_YUV444;
    }
    vlc_assert_unreachable();
}

static int CUDAAPI HandleVideoSequence(void *p_opaque, CUVIDEOFORMAT *p_format)
{
    decoder_t *p_dec = (decoder_t *) p_opaque;
    nvdec_ctx_t *p_ctx = p_dec->p_sys;
    int ret;

    // update vlc's output format using NVDEC parser's output
    ret = UpdateOutputFormat(p_dec, p_format);
    if (ret != VLC_SUCCESS)
        return 0;

    CUVIDDECODECREATEINFO dparams = {
        .ulWidth             = p_format->coded_width,
        .ulHeight            = p_format->coded_height,
        .ulTargetWidth       = p_dec->fmt_out.video.i_width,
        .ulTargetHeight      = p_dec->fmt_out.video.i_height,
        .bitDepthMinus8      = p_format->bit_depth_luma_minus8,
        .OutputFormat        = MapSurfaceFmt(p_dec->fmt_out.video.i_chroma),
        .CodecType           = p_format->codec,
        .ChromaFormat        = p_format->chroma_format,
        .ulNumDecodeSurfaces = p_ctx->i_nb_surface,
        .ulNumOutputSurfaces = 1,
        .DeinterlaceMode     = p_ctx->deintMode
    };

    ret = CUDA_CALL(p_ctx->cudaFunctions->cuCtxPushCurrent(p_ctx->cuCtx));
    if (ret != VLC_SUCCESS)
        return 0;

    /*if (p_ctx->cudecoder)
        CUDA_CALL(p_ctx->functions->cuvidDestroyDecoder(p_ctx->cudecoder));*/

    ret = CUDA_CALL(p_ctx->functions->cuvidCreateDecoder(&p_ctx->cudecoder, &dparams));
    CUDA_CALL(p_ctx->cudaFunctions->cuCtxPopCurrent(NULL));

    return (ret == VLC_SUCCESS);
}

static int CUDAAPI HandlePictureDecode(void *p_opaque, CUVIDPICPARAMS *p_picparams)
{
    decoder_t *p_dec = (decoder_t *) p_opaque;
    nvdec_ctx_t *p_ctx = p_dec->p_sys;
    int ret;

    ret = CUDA_CALL(p_ctx->cudaFunctions->cuCtxPushCurrent(p_ctx->cuCtx));
    if (ret != VLC_SUCCESS)
        return 0;

    ret = CUDA_CALL(p_ctx->functions->cuvidDecodePicture(p_ctx->cudecoder, p_picparams));
    CUDA_CALL(p_ctx->cudaFunctions->cuCtxPopCurrent(NULL));

    return (ret == VLC_SUCCESS);
}

static int CUDAAPI HandlePictureDisplay(void *p_opaque, CUVIDPARSERDISPINFO *p_dispinfo)
{
    decoder_t *p_dec = (decoder_t *) p_opaque;
    nvdec_ctx_t *p_ctx = p_dec->p_sys;
    int result;

    result = CUDA_CALL(p_ctx->cudaFunctions->cuCtxPushCurrent(p_ctx->cuCtx));
    if (result != VLC_SUCCESS)
        return 0;

    picture_t * p_pic = decoder_NewPicture(p_dec);

    CUdeviceptr cu_frame = 0;
    unsigned int i_pitch;
    CUVIDPROCPARAMS params;
    memset(&params, 0, sizeof(params));
    params.progressive_frame = p_ctx->b_progressive;
    params.top_field_first = p_dispinfo->top_field_first;

    // Map decoded frame to a device pointer
    result = CUDA_CALL(p_ctx->functions->cuvidMapVideoFrame(p_ctx->cudecoder, p_dispinfo->picture_index,
                                                            &cu_frame, &i_pitch, &params));
    if (result != VLC_SUCCESS)
        goto error;

    // Copy decoded frame into a new VLC picture
    p_pic->b_progressive = params.progressive_frame;
    p_pic->b_top_field_first = params.top_field_first;
    p_pic->date = p_dispinfo->timestamp;

    int i_plane_offset = 0, height = 0;
    for (int i_plane = 0; i_plane < p_pic->i_planes; i_plane++) {
        plane_t plane = p_pic->p[i_plane];
        height = p_dec->fmt_out.video.i_height >> ((i_plane) ? 1 : 0);
        CUDA_MEMCPY2D cu_cpy = {
            .srcMemoryType  = CU_MEMORYTYPE_DEVICE,
            .dstMemoryType  = CU_MEMORYTYPE_HOST,
            .srcDevice      = cu_frame,
            .srcY           = i_plane_offset,
            .srcPitch       = i_pitch,
            .dstHost        = plane.p_pixels,
            .dstPitch       = plane.i_pitch,
            .WidthInBytes   = __MIN(i_pitch, (size_t) plane.i_pitch),
            .Height         = height,
        };
        result = CUDA_CALL(p_ctx->cudaFunctions->cuMemcpy2D(&cu_cpy));
        if (result != VLC_SUCCESS)
            goto error;
        i_plane_offset += cu_cpy.Height;
    }

    // Release surface on GPU
    result = CUDA_CALL(p_ctx->functions->cuvidUnmapVideoFrame(p_ctx->cudecoder, cu_frame));
    if (result != VLC_SUCCESS)
        goto error;

    // Push decoded frame to display queue
    decoder_QueueVideo(p_dec, p_pic);

    CUDA_CALL(p_ctx->cudaFunctions->cuCtxPopCurrent(NULL));
    return 1;

error:
    CUDA_CALL(p_ctx->cudaFunctions->cuCtxPopCurrent(NULL));
    picture_Release(p_pic);
    return 0;
}

static int CuvidPushBlock(decoder_t *p_dec, block_t *p_block)
{
    nvdec_ctx_t *p_ctx = p_dec->p_sys;

    CUVIDSOURCEDATAPACKET cupacket = {0};
    cupacket.flags |= CUVID_PKT_TIMESTAMP;
    cupacket.payload_size = p_block->i_buffer;
    cupacket.payload = p_block->p_buffer;
    cupacket.timestamp = p_block->i_pts == VLC_TICK_INVALID ? p_block->i_dts : p_block->i_pts;

    return CUDA_CALL(p_ctx->functions->cuvidParseVideoData(p_ctx->cuparser, &cupacket));
}

static block_t * HXXXProcessBlock(decoder_t *p_dec, block_t *p_block)
{
    nvdec_ctx_t *p_ctx = p_dec->p_sys;
    if (p_ctx->hh.b_is_xvcC && !p_ctx->b_xps_pushed) {
        block_t *p_xps_blocks;   // parameter set blocks (SPS/PPS/VPS)
        if (p_dec->fmt_in.i_codec == VLC_CODEC_H264) {
            p_xps_blocks = h264_helper_get_annexb_config(&p_ctx->hh);
        } else if (p_dec->fmt_in.i_codec == VLC_CODEC_HEVC) {
            p_xps_blocks = hevc_helper_get_annexb_config(&p_ctx->hh);
        } else {
            return NULL;
        }
        for (block_t *p_b = p_xps_blocks; p_b != NULL; p_b = p_b->p_next) {
            CuvidPushBlock(p_dec, p_b);
        }
        p_ctx->b_xps_pushed = true;
    }

    return p_ctx->hh.pf_process_block(&p_ctx->hh, p_block, NULL);
}

static int CuvidPushEOS(decoder_t *p_dec)
{
    nvdec_ctx_t *p_ctx = p_dec->p_sys;

    CUVIDSOURCEDATAPACKET cupacket = {0};
    cupacket.flags |= CUVID_PKT_ENDOFSTREAM;
    cupacket.payload_size = 0;
    cupacket.payload = NULL;
    cupacket.timestamp = 0;

    return CUDA_CALL(p_ctx->functions->cuvidParseVideoData(p_ctx->cuparser, &cupacket));
}

static int DecodeBlock(decoder_t *p_dec, block_t *p_block)
{
    nvdec_ctx_t *p_ctx = p_dec->p_sys;
    if (p_block == NULL) {
        // Flush stream
        return CuvidPushEOS(p_dec);
    }
    if (p_ctx->b_is_hxxx) {
        p_block = HXXXProcessBlock(p_dec, p_block);
        if (p_block == NULL) {
            // try next block
            return VLCDEC_SUCCESS;
        }
    }
    if ((p_block->i_flags & BLOCK_FLAG_INTERLACED_MASK) && p_ctx->b_progressive) {
        p_ctx->b_progressive = false;
        p_ctx->deintMode = cudaVideoDeinterlaceMode_Adaptive;
    }
    return CuvidPushBlock(p_dec, p_block);
}

static int MapCodecID(int i_vlc_fourcc)
{
    switch (i_vlc_fourcc) {
        case VLC_CODEC_H264: return cudaVideoCodec_H264;
        case VLC_CODEC_HEVC: return cudaVideoCodec_HEVC;
        case VLC_CODEC_VP8:  return cudaVideoCodec_VP8;
        case VLC_CODEC_VP9:  return cudaVideoCodec_VP9;
    }
    vlc_assert_unreachable();
}

static int OpenDecoder(vlc_object_t *p_this)
{
    decoder_t *p_dec = (decoder_t *) p_this;
    nvdec_ctx_t *p_ctx;
    int result;
    p_ctx = calloc(1, sizeof(*p_ctx));
    if (!p_ctx)
        return VLC_ENOMEM;

    p_dec->p_sys = p_ctx;
    p_dec->pf_decode = DecodeBlock;

    switch (p_dec->fmt_in.i_codec) {
        case VLC_CODEC_H264:
        case VLC_CODEC_HEVC:
            p_ctx->b_is_hxxx = true;
            p_ctx->i_nb_surface = MAX_HXXX_SURFACES;
            hxxx_helper_init(&p_ctx->hh, VLC_OBJECT(p_dec),
                             p_dec->fmt_in.i_codec, false);
            result = hxxx_helper_set_extra(&p_ctx->hh, p_dec->fmt_in.p_extra,
                                           p_dec->fmt_in.i_extra);
            if (result != VLC_SUCCESS) {
                hxxx_helper_clean(&p_ctx->hh);
                free(p_ctx);
                return VLC_EGENERIC;
            }
            break;
        // VP8/VP9 are unsupported for now
        case VLC_CODEC_VP8:
        case VLC_CODEC_VP9:
        default:
            free(p_ctx);
            return VLC_EGENERIC;
    }

    result = cuvid_load_functions(&p_ctx->functions, NULL);
    if (result != VLC_SUCCESS) {
        hxxx_helper_clean(&p_ctx->hh);
        free(p_ctx);
        return VLC_EGENERIC;
    }
    result = cuda_load_functions(&p_ctx->cudaFunctions, NULL);
    if (result != VLC_SUCCESS) {
        hxxx_helper_clean(&p_ctx->hh);
        free(p_ctx);
        return VLC_EGENERIC;
    }

    result = CUDA_CALL(p_ctx->cudaFunctions->cuInit(0));
    if (result != VLC_SUCCESS)
        goto error;
    result = CUDA_CALL(p_ctx->cudaFunctions->cuCtxCreate(&p_ctx->cuCtx, 0, 0));
    if (result != VLC_SUCCESS)
        goto error;

    CUVIDDECODECAPS caps = {
        .eCodecType         = MapCodecID(p_dec->fmt_in.i_codec),
        .eChromaFormat      = cudaVideoChromaFormat_420,
        .nBitDepthMinus8    = 0    // FIXME: hardcoded to 8-bit for now
    };

    result = CUDA_CALL(p_ctx->cudaFunctions->cuCtxPushCurrent(p_ctx->cuCtx));
    if (result != VLC_SUCCESS)
        goto error;
    result =  CUDA_CALL(p_ctx->functions->cuvidGetDecoderCaps(&caps));
    CUDA_CALL(p_ctx->cudaFunctions->cuCtxPopCurrent(NULL));
    if (result != VLC_SUCCESS || !caps.bIsSupported) {
        msg_Err(p_dec, "No hardware for NVDEC");
        goto error;
    }

    p_ctx->b_progressive    = true;
    p_ctx->deintMode        = cudaVideoDeinterlaceMode_Weave;

    p_ctx->pparams.CodecType               = MapCodecID(p_dec->fmt_in.i_codec);
    p_ctx->pparams.ulClockRate             = CLOCK_FREQ;
    p_ctx->pparams.ulMaxDisplayDelay       = NVDEC_PIPELINE_SIZE;
    p_ctx->pparams.ulMaxNumDecodeSurfaces  = p_ctx->i_nb_surface;
    p_ctx->pparams.pUserData               = p_dec;
    p_ctx->pparams.pfnSequenceCallback     = HandleVideoSequence;
    p_ctx->pparams.pfnDecodePicture        = HandlePictureDecode;
    p_ctx->pparams.pfnDisplayPicture       = HandlePictureDisplay;
    result = CUDA_CALL(p_ctx->functions->cuvidCreateVideoParser(&p_ctx->cuparser, &p_ctx->pparams));
    if (result != VLC_SUCCESS) {
        msg_Err(p_dec, "Unable to create NVDEC video parser");
        goto error;
    }

    return VLC_SUCCESS;

error:
    CloseDecoder(p_this);
    return VLC_EGENERIC;
}

static void CloseDecoder(vlc_object_t *p_this)
{
    decoder_t *p_dec = (decoder_t *) p_this;
    nvdec_ctx_t *p_ctx = p_dec->p_sys;
    if (p_ctx->cudecoder)
        CUDA_CALL(p_ctx->functions->cuvidDestroyDecoder(p_ctx->cudecoder));
    if (p_ctx->cuparser)
        CUDA_CALL(p_ctx->functions->cuvidDestroyVideoParser(p_ctx->cuparser));
    if (p_ctx->cuCtx)
        CUDA_CALL(p_ctx->cudaFunctions->cuCtxDestroy(p_ctx->cuCtx));
    cuda_free_functions(&p_ctx->cudaFunctions);
    cuvid_free_functions(&p_ctx->functions);
    if (p_ctx->b_is_hxxx)
        hxxx_helper_clean(&p_ctx->hh);
    free(p_ctx);
}

