#include <../../../../../../AppData/Local/Programs/Python/Python35-32/include/Python.h>
#include <QQmlApplicationEngine>
#include <QApplication>
#include <QtQuick>
#include <QtCharts/QChartView>
#include <cmath>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

     QQuickView viewer;
     // The following are needed to make examples run without having to install the module
     // in desktop environments.
    #ifdef Q_OS_WIN
        QString extraImportPath(QStringLiteral("%1/../../../../%2"));
    #else
        QString extraImportPath(QStringLiteral("%1/../../../%2"));
    #endif

    //import path define in viewer ?
    viewer.engine()->addImportPath(extraImportPath.arg(QGuiApplication::applicationDirPath(),
                                       QString::fromLatin1("qml")));

    QObject::connect(viewer.engine(), &QQmlEngine::quit, &viewer, &QWindow::close);
    viewer.setTitle(QStringLiteral("QML Curve Input"));

    viewer.setSource(QUrl("qrc:/main.qml"));

    //Py_SetProgramName(L'ab');  /* optional but recommended */
    Py_Initialize();

    ///////////////////////////////////
    PyRun_SimpleString("result = 5 ** 2");
    PyObject * module = PyImport_AddModule("__main__"); // borrowed reference
    assert(module);                                     // __main__ should always exist
    PyObject * dictionary = PyModule_GetDict(module);   // borrowed reference
    assert(dictionary);                                 // __main__ should have a dictionary
    PyObject * result  = PyDict_GetItemString(dictionary, "result");     // borrowed reference

    assert(result);                                     // just added result
    long result_value = PyLong_AsLong(result);       // already checked that it is an int

    qDebug() << result_value;



    //////////////////////////////////

    Py_Finalize();


    return app.exec();
}

