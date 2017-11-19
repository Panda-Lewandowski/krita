#ifndef TESTCOLORBLENDING_H
#define TESTCOLORBLENDING_H

#include <QtTest>
#include "ui_TestColorBlendingWindow.h"
#include <cmath>
#include <limits>

#include "kritapigment_export.h"
#include <KoIntegerMaths.h>
#include "KoChannelInfo.h"
#include "KoLut.h"
#include "KoColorSpaceMaths.h"

class TestColorBlending : public QObject
{
    Q_OBJECT
private Q_SLOTS:
    void main_test_window();
};

class TestColorBlendingWindow: public QDialog, Ui::MainWindow
{
    Q_OBJECT
public:
    TestColorBlendingWindow(QWidget *parent = 0) : QDialog(parent) {
        setupUi(this);
    }
};


template < typename _T, typename _Tdst = _T >
class KoColorSpaceOptimizedMaths : public KoColorSpaceMaths<_T, _T>
{
    typedef KoColorSpaceMathsTraits<_T> traits;
    typedef typename traits::compositetype src_compositetype;
    typedef typename KoColorSpaceMathsTraits<_Tdst>::compositetype dst_compositetype;
public:
    inline static _Tdst multiply(_T a, _Tdst b);
    inline static _Tdst multiply(_T a, _Tdst b, _Tdst c);
    inline static dst_compositetype divide(_T a, _Tdst b);
    inline static _T invert(_T a);
    inline static _T blend(_T a, _T b, _T alpha);
    inline static _Tdst scaleToA(_T a);
    inline static dst_compositetype clamp(dst_compositetype val);
    inline static _Tdst clampAfterScale(dst_compositetype val);
};

#endif // TESTCOLORBLENDING_H
