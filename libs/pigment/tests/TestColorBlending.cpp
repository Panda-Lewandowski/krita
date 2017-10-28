#include "TestColorBlending.h"

#include <QTest>
#include <QDialog>


void TestColorBlending::test()
{
    QDialog window;
    window.exec();
}

TestColorBlendingWindow::TestColorBlendingWindow(QWidget *parent)
    : QWidget(parent)
{
   //this->setupUi(this);
}

QTEST_MAIN(TestColorBlending)
