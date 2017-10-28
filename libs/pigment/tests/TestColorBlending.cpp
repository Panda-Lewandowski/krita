#include "TestColorBlending.h"

#include <QTest>
#include <QApplication>


void TestColorBlending::test()
{
    //TestColorBlendingWindow w;
    //w.show();
    //QApplication app();

    TestColorBlendingWindow *demo = new TestColorBlendingWindow;
    demo->show();
    //demo->exec();
}



QTEST_MAIN(TestColorBlending)
/*int main(int argc, char **argv)
{
    QApplication app(argc, argv);

    TestColorBlendingWindow *demo = new TestColorBlendingWindow;
    demo->show();
    return app.exec();
}*/
