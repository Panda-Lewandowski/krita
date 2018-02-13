#include "TestColorBlending.h"

#include <QTest>
#include <KoColorSpace.h>
#include <KoCompositeOp.h>
#include <KoColorSpaceRegistry.h>
#include <KoColorSpaceTraits.h>
#include <kis_debug.h>

#if defined _MSC_VER
#define MEMALIGN_ALLOC(p, a, s) ((*(p)) = _aligned_malloc((s), (a)), *(p) ? 0 : errno)
#define MEMALIGN_FREE(p) _aligned_free((p))
#else
#define MEMALIGN_ALLOC(p, a, s) posix_memalign((p), (a), (s))
#define MEMALIGN_FREE(p) free((p))
#endif

const int alpha_pos = 3;

enum AlphaRange {
    ALPHA_ZERO,
    ALPHA_UNIT,
    ALPHA_RANDOM
};


template <typename channel_type, class RandomGenerator>
inline channel_type generateAlphaValue(AlphaRange range, RandomGenerator &rnd) {
    channel_type value = 0;

    switch (range) {
    case ALPHA_ZERO:
        break;
    case ALPHA_UNIT:
        value = rnd.unit();
        break;
    case ALPHA_RANDOM:
        value = rnd();
        break;
    }

    return value;
}

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/uniform_real.hpp>

template <typename channel_type>
struct RandomGenerator {
    channel_type operator() () {
        qFatal("Wrong template instantiation");
        return channel_type(0);
    }

    channel_type unit() {
        qFatal("Wrong template instantiation");
        return channel_type(0);
    }
};

template <>
struct RandomGenerator<quint8>
{
    RandomGenerator(int seed)
        : m_smallint(0,255),
          m_rnd(seed)
    {
    }

    quint8 operator() () {
        return m_smallint(m_rnd);
    }

    quint8 unit() {
        return KoColorSpaceMathsTraits<quint8>::unitValue;
    }

    boost::uniform_smallint<int> m_smallint;
    boost::mt11213b m_rnd;
};

template <>
struct RandomGenerator<float>
{
    RandomGenerator(int seed)
        : m_rnd(seed)
    {
    }

    float operator() () {
        //return float(m_rnd()) / float(m_rnd.max());
        return m_smallfloat(m_rnd);
    }

    float unit() {
        return KoColorSpaceMathsTraits<float>::unitValue;
    }

    boost::uniform_real<float> m_smallfloat;
    boost::mt11213b m_rnd;
};

template <>
struct RandomGenerator<double> : RandomGenerator<float>
{
    RandomGenerator(int seed)
        : RandomGenerator<float>(seed)
    {
    }
};


template <typename channel_type>
void generateDataLine(uint seed, int numPixels, quint8 *srcPixels, quint8 *dstPixels, quint8 *mask, AlphaRange srcAlphaRange, AlphaRange dstAlphaRange)
{
    Q_ASSERT(numPixels >= 4);

    RandomGenerator<channel_type> rnd(seed);
    RandomGenerator<quint8> maskRnd(seed + 1);

    channel_type *srcArray = reinterpret_cast<channel_type*>(srcPixels);
    channel_type *dstArray = reinterpret_cast<channel_type*>(dstPixels);

    for (int i = 0; i < numPixels; i++) {
        for (int j = 0; j < 3; j++) {
            channel_type s = rnd();
            channel_type d = rnd();
            *(srcArray++) = s;
            *(dstArray++) = d;
        }

        channel_type sa = generateAlphaValue<channel_type>(srcAlphaRange, rnd);
        channel_type da = generateAlphaValue<channel_type>(dstAlphaRange, rnd);
        *(srcArray++) = sa;
        *(dstArray++) = da;

        *(mask++) = maskRnd();
    }
}

void printData(int numPixels, quint8 *srcPixels, quint8 *dstPixels, quint8 *mask)
{
    for (int i = 0; i < numPixels; i++) {
        qDebug() << "Src: "
                 << srcPixels[i*4] << "\t"
                 << srcPixels[i*4+1] << "\t"
                 << srcPixels[i*4+2] << "\t"
                 << srcPixels[i*4+3] << "\t"
                 << "Msk:" << mask[i];

        qDebug() << "Dst: "
                 << dstPixels[i*4] << "\t"
                 << dstPixels[i*4+1] << "\t"
                 << dstPixels[i*4+2] << "\t"
                 << dstPixels[i*4+3];
    }
}

const int rowStride = 64;
const int totalRows = 64;
const QRect processRect(0,0,64,64);
const int numPixels = rowStride * totalRows;
const int numTiles = 1024;


struct Tile {
    quint8 *src;
    quint8 *dst;
    quint8 *mask;
};
#include <stdint.h>
QVector<Tile> generateTiles(int size,
                            const int srcAlignmentShift,
                            const int dstAlignmentShift,
                            AlphaRange srcAlphaRange,
                            AlphaRange dstAlphaRange,
                            const quint32 pixelSize)
{
    QVector<Tile> tiles(size);

#ifdef HAVE_VC
    const int vecSize = Vc::float_v::size();
#else
    const int vecSize = 1;
#endif

    // the 256 are used to make sure that we have a good alignment no matter what build options are used.
    const size_t pixelAlignment = qMax(size_t(vecSize * sizeof(float)), size_t(256));
    const size_t maskAlignment = qMax(size_t(vecSize), size_t(256));
    for (int i = 0; i < size; i++) {
        void *ptr = 0;
        int error = MEMALIGN_ALLOC(&ptr, pixelAlignment, numPixels * pixelSize + srcAlignmentShift);
        if (error) {
            qFatal("posix_memalign failed: %d", error);
        }
        tiles[i].src = (quint8*)ptr + srcAlignmentShift;
        error = MEMALIGN_ALLOC(&ptr, pixelAlignment, numPixels * pixelSize + dstAlignmentShift);
        if (error) {
            qFatal("posix_memalign failed: %d", error);
        }
        tiles[i].dst = (quint8*)ptr + dstAlignmentShift;
        error = MEMALIGN_ALLOC(&ptr, maskAlignment, numPixels);
        if (error) {
            qFatal("posix_memalign failed: %d", error);
        }
        tiles[i].mask = (quint8*)ptr;

        if (pixelSize == 4) {
            generateDataLine<quint8>(1, numPixels, tiles[i].src, tiles[i].dst, tiles[i].mask, srcAlphaRange, dstAlphaRange);
        } else if (pixelSize == 16) {
            generateDataLine<float>(1, numPixels, tiles[i].src, tiles[i].dst, tiles[i].mask, srcAlphaRange, dstAlphaRange);
        } else {
            qFatal("Pixel size %i is not implemented", pixelSize);
        }
    }

    return tiles;
}

void freeTiles(QVector<Tile> tiles,
               const int srcAlignmentShift,
               const int dstAlignmentShift)
{
    Q_FOREACH (const Tile &tile, tiles) {
        MEMALIGN_FREE(tile.src - srcAlignmentShift);
        MEMALIGN_FREE(tile.dst - dstAlignmentShift);
        MEMALIGN_FREE(tile.mask);
    }
}

template <typename channel_type>
inline bool fuzzyCompare(channel_type a, channel_type b, channel_type prec) {
    return qAbs(a - b) <= prec;
}

template <typename channel_type>
inline bool comparePixels(channel_type *p1, channel_type *p2, channel_type prec) {
    return (p1[3] == p2[3] && p1[3] == 0) ||
        (fuzzyCompare(p1[0], p2[0], prec) &&
         fuzzyCompare(p1[1], p2[1], prec) &&
         fuzzyCompare(p1[2], p2[2], prec) &&
         fuzzyCompare(p1[3], p2[3], prec));
}

template <typename channel_type>
bool compareTwoOpsPixels(QVector<Tile> &tiles, channel_type prec) {
    channel_type *dst1 = reinterpret_cast<channel_type*>(tiles[0].dst);
    channel_type *dst2 = reinterpret_cast<channel_type*>(tiles[1].dst);

    channel_type *src1 = reinterpret_cast<channel_type*>(tiles[0].src);
    channel_type *src2 = reinterpret_cast<channel_type*>(tiles[1].src);

    for (int i = 0; i < numPixels; i++) {
        if (!comparePixels<channel_type>(dst1, dst2, prec)) {
            qDebug() << "Wrong result:" << i;
            qDebug() << "Act: " << dst1[0] << dst1[1] << dst1[2] << dst1[3];
            qDebug() << "Exp: " << dst2[0] << dst2[1] << dst2[2] << dst2[3];
            qDebug() << "Dif: " << dst1[0] - dst2[0] << dst1[1] - dst2[1] << dst1[2] - dst2[2] << dst1[3] - dst2[3];

            channel_type *s1 = src1 + 4 * i;
            channel_type *s2 = src2 + 4 * i;

            qDebug() << "SrcA:" << s1[0] << s1[1] << s1[2] << s1[3];
            qDebug() << "SrcE:" << s2[0] << s2[1] << s2[2] << s2[3];

            qDebug() << "MskA:" << tiles[0].mask[i];
            qDebug() << "MskE:" << tiles[1].mask[i];

            return false;
        }
        dst1 += 4;
        dst2 += 4;
    }
    return true;
}

bool compareTwoOps(bool haveMask, const KoCompositeOp *op1, const KoCompositeOp *op2)
{
    Q_ASSERT(op1->colorSpace()->pixelSize() == op2->colorSpace()->pixelSize());
    const quint32 pixelSize = op1->colorSpace()->pixelSize();
    const int alignment = 16;
    QVector<Tile> tiles = generateTiles(2, alignment, alignment, ALPHA_RANDOM, ALPHA_RANDOM, op1->colorSpace()->pixelSize());

    KoCompositeOp::ParameterInfo params;
    params.dstRowStride  = 4 * rowStride;
    params.srcRowStride  = 4 * rowStride;
    params.maskRowStride = rowStride;
    params.rows          = processRect.height();
    params.cols          = processRect.width();
    // This is a hack as in the old version we get a rounding of opacity to this value
    params.opacity       = float(Arithmetic::scale<quint8>(0.5*1.0f))/255.0;
    params.flow          = 0.3*1.0f;
    params.channelFlags  = QBitArray();

    params.dstRowStart   = tiles[0].dst;
    params.srcRowStart   = tiles[0].src;
    params.maskRowStart  = haveMask ? tiles[0].mask : 0;
    op1->composite(params);

    params.dstRowStart   = tiles[1].dst;
    params.srcRowStart   = tiles[1].src;
    params.maskRowStart  = haveMask ? tiles[1].mask : 0;
    op2->composite(params);

    bool compareResult = true;
    if (pixelSize == 4) {
        compareResult = compareTwoOpsPixels<quint8>(tiles, 10);
    }
    else if (pixelSize == 16) {
        compareResult = compareTwoOpsPixels<float>(tiles, 2e-7);
    }
    else {
        qFatal("Pixel size %i is not implemented", pixelSize);
    }

    freeTiles(tiles, alignment, alignment);

    return compareResult;
}

QString getTestName(bool haveMask,
                    const int srcAlignmentShift,
                    const int dstAlignmentShift,
                    AlphaRange srcAlphaRange,
                    AlphaRange dstAlphaRange)
{

    QString testName;
    testName +=
        !srcAlignmentShift && !dstAlignmentShift ? "Aligned   " :
        !srcAlignmentShift &&  dstAlignmentShift ? "SrcUnalig " :
         srcAlignmentShift && !dstAlignmentShift ? "DstUnalig " :
         srcAlignmentShift &&  dstAlignmentShift ? "Unaligned " : "###";

    testName += haveMask ? "Mask   " : "NoMask ";

    testName +=
        srcAlphaRange == ALPHA_RANDOM ? "SrcRand " :
        srcAlphaRange == ALPHA_ZERO   ? "SrcZero " :
        srcAlphaRange == ALPHA_UNIT   ? "SrcUnit " : "###";

    testName +=
        dstAlphaRange == ALPHA_RANDOM ? "DstRand" :
        dstAlphaRange == ALPHA_ZERO   ? "DstZero" :
        dstAlphaRange == ALPHA_UNIT   ? "DstUnit" : "###";

    return testName;
}

void benchmarkCompositeOp(const KoCompositeOp *op,
                          bool haveMask,
                          qreal opacity,
                          qreal flow,
                          const int srcAlignmentShift,
                          const int dstAlignmentShift,
                          AlphaRange srcAlphaRange,
                          AlphaRange dstAlphaRange)
{
    QString testName = getTestName(haveMask, srcAlignmentShift, dstAlignmentShift, srcAlphaRange, dstAlphaRange);

    QVector<Tile> tiles =
        generateTiles(numTiles, srcAlignmentShift, dstAlignmentShift, srcAlphaRange, dstAlphaRange, op->colorSpace()->pixelSize());

    const int tileOffset = 4 * (processRect.y() * rowStride + processRect.x());

    KoCompositeOp::ParameterInfo params;
    params.dstRowStride  = 4 * rowStride;
    params.srcRowStride  = 4 * rowStride;
    params.maskRowStride = rowStride;
    params.rows          = processRect.height();
    params.cols          = processRect.width();
    params.opacity       = opacity;
    params.flow          = flow;
    params.channelFlags  = QBitArray();

    QTime timer;
    timer.start();

    Q_FOREACH (const Tile &tile, tiles) {
        params.dstRowStart   = tile.dst + tileOffset;
        params.srcRowStart   = tile.src + tileOffset;
        params.maskRowStart  = haveMask ? tile.mask : 0;
        op->composite(params);
    }

    qDebug() << testName << "RESULT:" << timer.elapsed() << "msec";

    freeTiles(tiles, srcAlignmentShift, dstAlignmentShift);
}

void benchmarkCompositeOp(const KoCompositeOp *op, const QString &postfix)
{
    qDebug() << "Testing Composite Op:" << op->id() << "(" << postfix << ")";

    benchmarkCompositeOp(op, true, 0.5, 0.3, 0, 0, ALPHA_RANDOM, ALPHA_RANDOM);
    benchmarkCompositeOp(op, true, 0.5, 0.3, 8, 0, ALPHA_RANDOM, ALPHA_RANDOM);
    benchmarkCompositeOp(op, true, 0.5, 0.3, 0, 8, ALPHA_RANDOM, ALPHA_RANDOM);
    benchmarkCompositeOp(op, true, 0.5, 0.3, 4, 8, ALPHA_RANDOM, ALPHA_RANDOM);

/// --- Vary the content of the source and destination

    benchmarkCompositeOp(op, false, 1.0, 1.0, 0, 0, ALPHA_RANDOM, ALPHA_RANDOM);
    benchmarkCompositeOp(op, false, 1.0, 1.0, 0, 0, ALPHA_ZERO, ALPHA_RANDOM);
    benchmarkCompositeOp(op, false, 1.0, 1.0, 0, 0, ALPHA_UNIT, ALPHA_RANDOM);

/// ---

    benchmarkCompositeOp(op, false, 1.0, 1.0, 0, 0, ALPHA_RANDOM, ALPHA_ZERO);
    benchmarkCompositeOp(op, false, 1.0, 1.0, 0, 0, ALPHA_ZERO, ALPHA_ZERO);
    benchmarkCompositeOp(op, false, 1.0, 1.0, 0, 0, ALPHA_UNIT, ALPHA_ZERO);

/// ---

    benchmarkCompositeOp(op, false, 1.0, 1.0, 0, 0, ALPHA_RANDOM, ALPHA_UNIT);
    benchmarkCompositeOp(op, false, 1.0, 1.0, 0, 0, ALPHA_ZERO, ALPHA_UNIT);
    benchmarkCompositeOp(op, false, 1.0, 1.0, 0, 0, ALPHA_UNIT, ALPHA_UNIT);
}



void TestColorBlending::init() {
    qApp->setApplicationName("krita");
    // All Krita's resource types
    KoResourcePaths::addResourceType("kis_pics", "data", "/pics/");
    KoResourcePaths::addResourceType("kis_images", "data", "/images/");
    KoResourcePaths::addResourceType("icc_profiles", "data", "/profiles/");
    KoResourcePaths::addResourceType("metadata_schema", "data", "/metadata/schemas/");
    KoResourcePaths::addResourceType("kis_brushes", "data", "/brushes/");
    KoResourcePaths::addResourceType("kis_taskset", "data", "/taskset/");
    KoResourcePaths::addResourceType("kis_taskset", "data", "/taskset/");
    KoResourcePaths::addResourceType("gmic_definitions", "data", "/gmic/");
    KoResourcePaths::addResourceType("kis_resourcebundles", "data", "/bundles/");
    KoResourcePaths::addResourceType("kis_defaultpresets", "data", "/defaultpresets/");
    KoResourcePaths::addResourceType("kis_paintoppresets", "data", "/paintoppresets/");
    KoResourcePaths::addResourceType("kis_workspaces", "data", "/workspaces/");
    KoResourcePaths::addResourceType("psd_layer_style_collections", "data", "/asl");
    KoResourcePaths::addResourceType("ko_patterns", "data", "/patterns/", true);
    KoResourcePaths::addResourceType("ko_gradients", "data", "/gradients/");
    KoResourcePaths::addResourceType("ko_gradients", "data", "/gradients/", true);
    KoResourcePaths::addResourceType("ko_palettes", "data", "/palettes/", true);
    KoResourcePaths::addResourceType("kis_shortcuts", "data", "/shortcuts/");
    KoResourcePaths::addResourceType("kis_actions", "data", "/actions");
    KoResourcePaths::addResourceType("icc_profiles", "data", "/color/icc");
    KoResourcePaths::addResourceType("ko_effects", "data", "/effects/");
    KoResourcePaths::addResourceType("tags", "data", "/tags/");
    KoResourcePaths::addResourceType("templates", "data", "/templates");
    KoResourcePaths::addResourceType("pythonscripts", "data", "/pykrita");
    KoResourcePaths::addResourceType("symbols", "data", "/symbols");

    //    // Extra directories to look for create resources. (Does anyone actually use that anymore?)
    //    KoResourcePaths::addResourceDir("ko_gradients", "/usr/share/create/gradients/gimp");
    //    KoResourcePaths::addResourceDir("ko_gradients", QDir::homePath() + QString("/.create/gradients/gimp"));
    //    KoResourcePaths::addResourceDir("ko_patterns", "/usr/share/create/patterns/gimp");
    //    KoResourcePaths::addResourceDir("ko_patterns", QDir::homePath() + QString("/.create/patterns/gimp"));
    //    KoResourcePaths::addResourceDir("kis_brushes", "/usr/share/create/brushes/gimp");
    //    KoResourcePaths::addResourceDir("kis_brushes", QDir::homePath() + QString("/.create/brushes/gimp"));
    //    KoResourcePaths::addResourceDir("ko_palettes", "/usr/share/create/swatches");
    //    KoResourcePaths::addResourceDir("ko_palettes", QDir::homePath() + QString("/.create/swatches"));

    // Make directories for all resources we can save, and tags
    QDir d;
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/tags/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/asl/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/bundles/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/gradients/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/paintoppresets/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/palettes/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/patterns/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/taskset/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/workspaces/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/input/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/pykrita/");
    d.mkpath(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/symbols/");

    // Indicate that it is now safe for users of KoResourcePaths to load resources
    KoResourcePaths::setReady();

}

//void TestColorBlending::benchmark()
//{
//    const KoColorSpace *cs = KoColorSpaceRegistry::instance()->rgb8();
//    KoCompositeOp *op = new KoDoubleOptimizedCompositeOpOver32<Vc::CurrentImplementation::current()>(cs);
//    benchmarkCompositeOp(op, "DoupleOptimized");
//    delete op;
//}

//void TestColorBlending::compareOverOps()
//{
//    const KoColorSpace *cs = KoColorSpaceRegistry::instance()->rgb8();
//    KoCompositeOp *opAct = KoOptimizedCompositeOpFactory::createOverOp32(cs);
//    KoCompositeOp *opExp = new KoDoubleOptimizedCompositeOpOver32<Vc::CurrentImplementation::current()>(cs);

//    QVERIFY(compareTwoOps(true, opAct, opExp));

//    delete opExp;
//    delete opAct;
//}

//void TestColorBlending::compareOverOpsNoMask()
//{
//    const KoColorSpace *cs = KoColorSpaceRegistry::instance()->rgb8();
//    KoCompositeOp *opAct = KoOptimizedCompositeOpFactory::createOverOp32(cs);
//    KoCompositeOp *opExp = new KoDoubleOptimizedCompositeOpOver32<Vc::CurrentImplementation::current()>(cs);

//    QVERIFY(compareTwoOps(false, opAct, opExp));

//    delete opExp;
//    delete opAct;
//}


void TestColorBlending::mask_test()
{
    Vc::uint8_v v;
    v[0] = 56; v[1] = 240; v[2] = 38; v[3] = 70;
    Vc::uint16_v mask1 = KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_mask_8_uint16((quint8*)&v);
    Vc::float_v mask2 = KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_mask_8((quint8*)&v);
    qDebug() << "Mymask:" << mask1[0] << mask1[1] << mask1[2] << mask1[3];
    qDebug() << "Optmask:" << mask2[0] << mask2[1] << mask2[2] << mask2[3];

}



void TestColorBlending::alpha_test()
{
    Vc::uint8_v v;
    v[0] = 56; v[1] = 240; v[2] = 38; v[3] = 80;
    Vc::uint16_v mask1 = KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_alpha_uint16<true>((quint8*)&v);
    Vc::float_v mask2 = KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_alpha_32<true>((quint8*)&v);
    qDebug() << "Myalpha:" << mask1[0] << mask1[1] << mask1[2] << mask1[3];
    qDebug() << "Optalpha:" << mask2[0] << mask2[1] << mask2[2] << mask2[3];

    mask1 = KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_alpha_uint16<false>((quint8*)&v);
    mask2 = KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_alpha_32<false>((quint8*)&v);
    qDebug() << "Myalpha:" << mask1[0] << mask1[1] << mask1[2] << mask1[3];
    qDebug() << "Optalpha:" << mask2[0] << mask2[1] << mask2[2] << mask2[3];

}


void TestColorBlending::colors_test()
{
    Vc::uint8_v v;
    v[0] = 56;
    v[1] = 240;
    v[2] = 38;
    v[3] = 80;
    Vc::float_v src_c1;
    Vc::float_v src_c2;
    Vc::float_v src_c3;

    Vc::uint16_v src_c4;
    Vc::uint16_v src_c5;
    Vc::uint16_v src_c6;
    KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_colors_32<true>((quint8*)&v, src_c1, src_c2, src_c3);
    qDebug() << "OptColor1:" << src_c1[0] << src_c1[1] << src_c1[2] << src_c1[3];
    qDebug() << "OptColor2:" << src_c2[0] << src_c2[1] << src_c2[2] << src_c2[3];
    qDebug() << "OptColor3:" << src_c3[0] << src_c3[1] << src_c3[2] << src_c3[3];

    KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_colors_uint16<true>((quint8*)&v, src_c4, src_c5, src_c6);
    qDebug() << "MyColor1:" << src_c4[0] << src_c4[1] << src_c4[2] << src_c4[3];
    qDebug() << "MyColor2:" << src_c5[0] << src_c5[1] << src_c5[2] << src_c5[3];
    qDebug() << "MyColor3:" << src_c6[0] << src_c6[1] << src_c6[2] << src_c6[3];

    KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_colors_32<false>((quint8*)&v, src_c1, src_c2, src_c3);
    qDebug() << "OptColor1:" << src_c1[0] << src_c1[1] << src_c1[2] << src_c1[3];
    qDebug() << "OptColor2:" << src_c2[0] << src_c2[1] << src_c2[2] << src_c2[3];
    qDebug() << "OptColor3:" << src_c3[0] << src_c3[1] << src_c3[2] << src_c3[3];

    KoStreamedMath<Vc::CurrentImplementation::current()>::fetch_colors_uint16<false>((quint8*)&v, src_c4, src_c5, src_c6);
    qDebug() << "MyColor1:" << src_c4[0] << src_c4[1] << src_c4[2] << src_c4[3];
    qDebug() << "MyColor2:" << src_c5[0] << src_c5[1] << src_c5[2] << src_c5[3];
    qDebug() << "MyColor3:" << src_c6[0] << src_c6[1] << src_c6[2] << src_c6[3];



}


//void TestColorBlending::channels_test()
//{
//    Vc::float_v src_alpha;src_alpha[0] = 46; src_alpha[1] = 0; src_alpha[2] = 0; src_alpha[3] = 0;
//    Vc::float_v src_c1;src_c1[0] = 54; src_c1[1] = 0; src_c1[2] = 0; src_c1[3] = 0;
//    Vc::float_v src_c2;src_c2[0] = 144; src_c2[1] = 0; src_c2[2] = 0; src_c2[3] = 0;
//    Vc::float_v src_c3;src_c3[0] = 233; src_c3[1] = 0; src_c3[2] = 0; src_c3[3] = 0;

//    Vc::uint16_v dst_alpha;dst_alpha[0] = 46; dst_alpha[1] = 0; dst_alpha[2] = 0; dst_alpha[3] = 0;
//    Vc::uint16_v dst_c1;dst_c1[0] = 54; dst_c1[1] = 0; dst_c1[2] = 0; dst_c1[3] = 0;
//    Vc::uint16_v dst_c2;dst_c2[0] = 144; dst_c2[1] = 0; dst_c2[2] = 0; dst_c2[3] = 0;
//    Vc::uint16_v dst_c3;dst_c3[0] = 233; dst_c3[1] = 0; dst_c3[2] = 0; dst_c3[3] = 0;

//    quint8* v;
//    quint8* w;

//    KoStreamedMath<Vc::CurrentImplementation::current()>::write_channels_32(v, src_alpha, src_c1, src_c2, src_c3);
//    qDebug() << "OptChannel:" << v[0] << v[1] << v[2] << v[3];

//    KoStreamedMath<Vc::CurrentImplementation::current()>::write_channels_uint16(w, dst_alpha,dst_c1, dst_c2, dst_c3);
//    //qDebug() << "MyChannel:" << w[0] << w[1] << w[2] << w[3];



//}

//void TestColorBlending::test_blend()
//{




//}

QTEST_MAIN(TestColorBlending)

