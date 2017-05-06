#include <QApplication>
#include <QQmlApplicationEngine>
#include <QtWidgets/QApplication>
#include <QtCore/QDir>
#include "mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    MainWindow w;
    w.show();

    return app.exec();
}
