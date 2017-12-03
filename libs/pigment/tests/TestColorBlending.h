#ifndef TESTCOLORBLENDING_H
#define TESTCOLORBLENDING_H

#include <QtTest>
#include "ui_TestColorBlendingWindow.h"
#include <cmath>
#include <limits>

#include "kritapigment_export.h"
#include "KoCompositeOpBase.h"
#include "KoCompositeOpRegistry.h"
#include "KoStreamedMath.h"
#include "KoOptimizedCompositeOpOver32.h"
#include <kis_debug.h>
#include "KoResourcePaths.h"


class TestColorBlending : public QObject
{
    Q_OBJECT
private Q_SLOTS:
    void init();

    //void main_test_window();
    void test();


};

class TestColorBlendingWindow: public QDialog, Ui::MainWindow
{
    Q_OBJECT
public:
    TestColorBlendingWindow(QWidget *parent = 0) : QDialog(parent) {
        setupUi(this);
    }
};

template<typename channels_type, typename pixel_type, bool alphaLocked, bool allChannelsFlag>
struct OptimizedOverCompositor32 {
    struct OptionalParams {
        OptionalParams(const KoCompositeOp::ParameterInfo& params)
            : channelFlags(params.channelFlags)
        {
        }
        const QBitArray &channelFlags;
    };

    template<bool haveMask, bool src_aligned, Vc::Implementation _impl>
    static ALWAYS_INLINE void compositeVector(const quint8 *src,
                                              quint8 *dst, const quint8 *mask,
                                              float opacity, const OptionalParams &oparams)
    {
        Q_UNUSED(oparams);

        Vc::float_v src_alpha;
        Vc::float_v dst_alpha;

        src_alpha = KoStreamedMath<_impl>::template fetch_alpha_32<src_aligned>(src);

        bool haveOpacity = opacity != 1.0;
        Vc::float_v opacity_norm_vec(opacity);

        Vc::float_v uint8Max((float)255.0);
        Vc::float_v uint8MaxRec1((float)1.0 / 255.0);
        Vc::float_v zeroValue(Vc::Zero);
        Vc::float_v oneValue(Vc::One);

        src_alpha *= opacity_norm_vec;

        if (haveMask) {
            Vc::float_v mask_vec = KoStreamedMath<_impl>::fetch_mask_8(mask);
            src_alpha *= mask_vec * uint8MaxRec1;
        }

        // The source cannot change the colors in the destination,
        // since its fully transparent
        if ((src_alpha == zeroValue).isFull()) {
            return;
        }

        dst_alpha = KoStreamedMath<_impl>::template fetch_alpha_32<true>(dst);

        Vc::float_v src_c1;
        Vc::float_v src_c2;
        Vc::float_v src_c3;

        Vc::float_v dst_c1;
        Vc::float_v dst_c2;
        Vc::float_v dst_c3;


        KoStreamedMath<_impl>::template fetch_colors_32<src_aligned>(src, src_c1, src_c2, src_c3);
        Vc::float_v src_blend;
        Vc::float_v new_alpha;

        if ((dst_alpha == uint8Max).isFull()) {
            new_alpha = dst_alpha;
            src_blend = src_alpha * uint8MaxRec1;
        } else if ((dst_alpha == zeroValue).isFull()) {
            new_alpha = src_alpha;
            src_blend = oneValue;
        } else {
            /**
             * The value of new_alpha can have *some* zero values,
             * which will result in NaN values while division. But
             * when converted to integers these NaN values will
             * be converted to zeroes, which is exactly what we need
             */
            new_alpha = dst_alpha + (uint8Max - dst_alpha) * src_alpha * uint8MaxRec1;

            // Optimized version of:
            //     src_blend = src_alpha / new_alpha;
            src_blend = OptiDiv<_impl>::divVector(src_alpha, new_alpha);

        }

        if (!(src_blend == oneValue).isFull()) {
            KoStreamedMath<_impl>::template fetch_colors_32<true>(dst, dst_c1, dst_c2, dst_c3);

            dst_c1 = src_blend * (src_c1 - dst_c1) + dst_c1;
            dst_c2 = src_blend * (src_c2 - dst_c2) + dst_c2;
            dst_c3 = src_blend * (src_c3 - dst_c3) + dst_c3;

        } else {
            if (!haveMask && !haveOpacity) {
                memcpy(dst, src, 4 * Vc::float_v::size());
                return;
            } else {
                // opacity has changed the alpha of the source,
                // so we can't just memcpy the bytes
                dst_c1 = src_c1;
                dst_c2 = src_c2;
                dst_c3 = src_c3;
            }
        }

        KoStreamedMath<_impl>::write_channels_32(dst, new_alpha, dst_c1, dst_c2, dst_c3);

    }

    template <bool haveMask, Vc::Implementation _impl>
    static ALWAYS_INLINE void compositeOnePixelScalar(const channels_type *src,
                                                      channels_type *dst, const quint8 *mask,
                                                      float opacity, const OptionalParams &oparams)
    {
        using namespace Arithmetic;
        const qint32 alpha_pos = 3;

        const float uint8Rec1 = 1.0 / 255.0;
        const float uint8Max = 255.0;

        float srcAlpha = src[alpha_pos];
        srcAlpha *= opacity;

        if (haveMask) {
            srcAlpha *= float(*mask) * uint8Rec1;
        }

        if (srcAlpha != 0.0) {

            float dstAlpha = dst[alpha_pos];
            float srcBlendNorm;

            if (dstAlpha == uint8Max) {
                srcBlendNorm = srcAlpha * uint8Rec1;
            } else if (dstAlpha == 0.0) {
                dstAlpha = srcAlpha;
                srcBlendNorm = 1.0;

                if (!allChannelsFlag) {
                    pixel_type *d = reinterpret_cast<pixel_type*>(dst);
                    *d = 0; // dstAlpha is already null
                }
            } else {
                dstAlpha += (uint8Max - dstAlpha) * srcAlpha * uint8Rec1;

                // Optimized version of:
                //     srcBlendNorm = srcAlpha / dstAlpha);
                srcBlendNorm = OptiDiv<_impl>::divScalar(srcAlpha, dstAlpha);


            }

            if(allChannelsFlag) {
                if (srcBlendNorm == 1.0) {
                    if (!alphaLocked) {
                        const pixel_type *s = reinterpret_cast<const pixel_type*>(src);
                        pixel_type *d = reinterpret_cast<pixel_type*>(dst);
                        *d = *s;
                    } else {
                        dst[0] = src[0];
                        dst[1] = src[1];
                        dst[2] = src[2];
                    }
                } else if (srcBlendNorm != 0.0){
                    dst[0] = KoStreamedMath<_impl>::lerp_mixed_u8_float(dst[0], src[0], srcBlendNorm);
                    dst[1] = KoStreamedMath<_impl>::lerp_mixed_u8_float(dst[1], src[1], srcBlendNorm);
                    dst[2] = KoStreamedMath<_impl>::lerp_mixed_u8_float(dst[2], src[2], srcBlendNorm);
                }
            } else {
                const QBitArray &channelFlags = oparams.channelFlags;

                if (srcBlendNorm == 1.0) {
                    if(channelFlags.at(0)) dst[0] = src[0];
                    if(channelFlags.at(1)) dst[1] = src[1];
                    if(channelFlags.at(2)) dst[2] = src[2];
                } else if (srcBlendNorm != 0.0) {
                    if(channelFlags.at(0)) dst[0] = KoStreamedMath<_impl>::lerp_mixed_u8_float(dst[0], src[0], srcBlendNorm);
                    if(channelFlags.at(1)) dst[1] = KoStreamedMath<_impl>::lerp_mixed_u8_float(dst[1], src[1], srcBlendNorm);
                    if(channelFlags.at(2)) dst[2] = KoStreamedMath<_impl>::lerp_mixed_u8_float(dst[2], src[2], srcBlendNorm);
                }
            }

            if (!alphaLocked) {
                dst[alpha_pos] = KoStreamedMath<_impl>::round_float_to_uint(dstAlpha);
            }
        }
    }

};

/**
 * An optimized version of a composite op for the use in 4 byte
 * colorspaces with alpha channel placed at the last byte of
 * the pixel: C1_C2_C3_A.
 */
template<Vc::Implementation _impl>
class KoDoubleOptimizedCompositeOpOver32 : public KoCompositeOp
{
public:
    KoDoubleOptimizedCompositeOpOver32(const KoColorSpace* cs)
        : KoCompositeOp(cs, COMPOSITE_OVER, i18n("Normal"), KoCompositeOp::categoryMix()) {}

    using KoCompositeOp::composite;

    virtual void composite(const KoCompositeOp::ParameterInfo& params) const
    {
        if(params.maskRowStart) {
            composite<true>(params);
        } else {
            composite<false>(params);
        }
    }

    template <bool haveMask>
    inline void composite(const KoCompositeOp::ParameterInfo& params) const {
        if (params.channelFlags.isEmpty() ||
            params.channelFlags == QBitArray(4, true)) {

            KoStreamedMath<_impl>::template genericComposite32<haveMask, false, OverCompositor32<quint8, quint32, false, true> >(params);
        } else {
            const bool allChannelsFlag =
                params.channelFlags.at(0) &&
                params.channelFlags.at(1) &&
                params.channelFlags.at(2);

            const bool alphaLocked =
                !params.channelFlags.at(3);

            if (allChannelsFlag && alphaLocked) {
                KoStreamedMath<_impl>::template genericComposite32_novector<haveMask, false, OptimizedOverCompositor32<quint8, quint32, true, true> >(params);
            } else if (!allChannelsFlag && !alphaLocked) {
                KoStreamedMath<_impl>::template genericComposite32_novector<haveMask, false, OptimizedOverCompositor32<quint8, quint32, false, false> >(params);
            } else /*if (!allChannelsFlag && alphaLocked) */{
                KoStreamedMath<_impl>::template genericComposite32_novector<haveMask, false, OptimizedOverCompositor32<quint8, quint32, true, false> >(params);
            }
        }
    }
};


#endif // TESTCOLORBLENDING_H
