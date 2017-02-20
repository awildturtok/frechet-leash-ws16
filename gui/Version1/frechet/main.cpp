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

    QQuickView viewer;

    // run Python script with system call
    /*std::string filename = "discreteFrechet.py";
    std::string command = "python ";
    command += filename;
    system(command.c_str());

    QProcess process;
    process.start("python", QStringList() << "C:/Users/mertjose/Documents/Repositories/Frechet_Distance/qml_chart_example/discreteFrechet.py");
    */

//    QObject *parent;
//    QString program = "cmd.exe";
//    QStringList arguments;
//    arguments << "/c python C:/Users/mertjose/Documents/Repositories/Frechet_Distance/qml_chart_example/discreteFrechet.py";

//    QProcess *myProcess = new QProcess(parent);
//    myProcess->start(program, arguments);

    // The following are needed to make examples run without having to install the module
    // in desktop environments.
//    #ifdef Q_OS_WIN
//        QString extraImportPath(QStringLiteral("%1/../../../../%2"));
//    #else
//        QString extraImportPath(QStringLiteral("%1/../../../%2"));
//    #endif



    //import path define in viewer ?
//    viewer.engine()->addImportPath(extraImportPath.arg(QGuiApplication::applicationDirPath(),
//                                      QString::fromLatin1("qml")));

    //QObject::connect(viewer.engine(), &QQmlEngine::quit, &viewer, &QWindow::close);
    //viewer.setTitle(QStringLiteral("QML Curve Input"));

    //QWidget container = QWidget::createWindowContainer(&viewer,this);
    //Testdata testData(&viewer);
    //viewer.rootContext()->setContextProperty("testData", &testData);
    //viewer.setSource(QUrl("qrc:/main.qml"));

    MainWindow w;
    w.show();

    return app.exec();
}
