#include "TestColorBlending.h"

#include <QTest>
#include <QDialog>


void TestColorBlending::test()
{
    QDialog window;

    TestColorBlendingWindow ui;
    ui.setupUi(window);

    window.exec();
}

QTEST_MAIN(TestColorBlending)
