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



QTEST_MAIN(TestColorBlending)

