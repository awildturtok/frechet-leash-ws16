#include <QApplication>
#include <QQmlApplicationEngine>
#include <QtQuick>
#include <QtCharts/QChartView>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QQmlApplicationEngine engine;
    engine.load(QUrl(QStringLiteral("qrc:/main.qml")));
    QObject *window = engine.rootObjects().at(0);
    QObject *chart = window->children().at(0);

    /*transform coordinates ?! http://www.qcustomplot.com/index.php/support/forum/93
     * http://doc.qt.io/qt-5/coordsys.html#window-viewport-conversion
     * !!! http://stackoverflow.com/questions/18140446/display-the-plot-values-on-mouse-over-detect-scatter-points
     *
     */
    //chart->

    return app.exec();
}

