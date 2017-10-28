#ifndef TESTCOLORBLENDING_H
#define TESTCOLORBLENDING_H

#include <QtTest>
#include "ui_TestColorBlendingWindow.h"

class TestColorBlending : public QObject
{
    Q_OBJECT
private Q_SLOTS:
    void test();
};

class TestColorBlendingWindow: public QWidget, Ui::MainWindow
{
    Q_OBJECT
public:
    TestColorBlendingWindow(QWidget* parent = 0);
};
#endif // TESTCOLORBLENDING_H
