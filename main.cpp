#include "MainSBI.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QString platformPluginPath = QCoreApplication::applicationDirPath() + "/platforms";
    QCoreApplication::addLibraryPath(platformPluginPath);

    MainSBI w;
    w.show();

    return a.exec();
}
