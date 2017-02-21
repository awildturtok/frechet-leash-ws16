#include <QApplication>
#include <QQmlApplicationEngine>
#include <QtQuick>
#include <QtCharts/QChartView>
#include <QtWidgets/QApplication>
#include <QtQml/QQmlContext>
#include <QtQuick/QQuickView>
#include <QtQml/QQmlEngine>
#include <QtCore/QDir>
#include "testdata.h"
#include "mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    qDebug() << app.applicationDirPath();

    QQuickView viewer;

    MainWindow w;
    w.show();

    return app.exec();
}
