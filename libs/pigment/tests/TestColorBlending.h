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

class TestColorBlendingWindow: public QMainWindow, Ui::MainWindow
{
    Q_OBJECT
public:
    TestColorBlendingWindow(QMainWindow *parent = 0) : QMainWindow(parent) {
        setupUi(this);
    }
};
#endif // TESTCOLORBLENDING_H
