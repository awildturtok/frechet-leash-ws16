
#include <QtCharts>
#include <QtQuick>
//#include <QtQuick/QQuickView>
//#include <QtQuick/QQuickItem>
#include <QtCore>
//#include <QtCore/QDebug>
//#include <QtCore/QtMath>
//#include <QtCore/QObject>
#include "testdata.h"
#include <iostream>
#include <fstream>

using namespace std;
string delimiter = ";";
ofstream myfile;

Testdata::Testdata(QQuickView *appViewer, QObject *parent):
    QObject(parent), m_appViewer(appViewer), m_index(-1)
{
    qRegisterMetaType<QAbstractSeries*>();
    qRegisterMetaType<QAbstractAxis*>();

    generateLineTyp(2, 5, 100);

}

void Testdata::update(QAbstractSeries *series)
{
    if (series) {
        QXYSeries *xySeries = static_cast<QXYSeries *>(series);
        m_index++;
        if (m_index > m_data.count() - 1)
            m_index = 0;

        QVector<QPointF> points = m_data.at(m_index);
        // Use replace instead of clear + append, it's optimized for performance
        xySeries->replace(points);
    }
}
void Testdata::cleanDataFile(){
    //clean data file first
    myfile.open("data.csv", std::ofstream::out | std::ofstream::trunc);
    myfile.close();
}

void Testdata::printPointSeries(QAbstractSeries *series)
{
    if (series) {
        QVector<QPointF> mypoints;
        QXYSeries *xySeries = static_cast<QXYSeries *>(series);
        string seriesId = xySeries->name().toStdString();
//        m_index++;
//        if (m_index > m_data.count() - 1)
//            m_index = 0;
        int pointcount=0;

        //myfile.open ("data.csv");
       //append to file
        myfile.open("data.csv", ios::out | ios::app );
        int seriescount=xySeries->count();
        myfile << "x"<< seriesId << delimiter << "y" << seriesId << endl;
        int i = 0;
        for(i = 0; i< xySeries->count();i++){
            QPointF currentpoint = xySeries->at(i);
            float x = currentpoint.x();
            float y = currentpoint.y();

            qDebug()<< "x: " << x << " y:" << y;
            myfile << x << delimiter << y << endl;

            pointcount++;
        }
         qDebug()<< "finish writing";
        myfile.close();
        //QVector<QPointF> points = m_data.at(0);
        // Use replace instead of clear + append, it's optimized for performance
       // xySeries->replace(points);
    }
}

void Testdata::generateLineTyp(int type, int rowCount, int colCount)
{
    foreach (QVector<QPointF> row, m_data)
        row.clear();
    m_data.clear();

    // Append the new data depending on the type
    for (int i(0); i < rowCount; i++) {
        QVector<QPointF> points;
        points.reserve(colCount);
        for (int j(0); j < colCount; j++) {
            qreal x(0);
            qreal y(0);
            switch (type) {
            case 0:
                // data with sin + random component
                y = qSin(3.14159265358979 / 50 * j) + 0.5 ;
                x = j;
                break;
            case 1:
                // linear data
                x = j;
                y = (qreal) i / 10;
                break;
            case 2:
                x=j;
                y=j*j;
            default:
                // unknown, do nothing
                break;
            }
            points.append(QPointF(x, y));
        }
        m_data.append(points);
}
}
