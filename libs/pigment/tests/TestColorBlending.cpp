#include "TestColorBlending.h"

#include <QTest>
#include <QApplication>


void TestColorBlending::main_test_window()
{
    QEventLoop loop;

    TestColorBlendingWindow *demo = new TestColorBlendingWindow;
    demo->show();
    loop.exec();
}


//------------------------------ quint8 specialization ------------------------------//

template<>
inline quint8 KoColorSpaceOptimizedMaths<quint8>::multiply(quint8 a, quint8 b)
{
    return (quint8)UINT8_MULT(a, b);
}


template<>
inline quint8 KoColorSpaceOptimizedMaths<quint8>::multiply(quint8 a, quint8 b, quint8 c)
{
    return (quint8)UINT8_MULT3(a, b, c);
}

template<>
inline KoColorSpaceMathsTraits<quint8>::compositetype
KoColorSpaceOptimizedMaths<quint8>::divide(quint8 a, quint8 b)
{
    return UINT8_DIVIDE(a, b);
}

template<>
inline quint8 KoColorSpaceOptimizedMaths<quint8>::invert(quint8 a)
{
    return ~a;
}

template<>
inline quint8 KoColorSpaceOptimizedMaths<quint8>::blend(quint8 a, quint8 b, quint8 c)
{
    return UINT8_BLEND(a, b, c);
}


QTEST_MAIN(TestColorBlending)

