#include "MainSBI.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainSBI w;
    w.show();

    return a.exec();
}
