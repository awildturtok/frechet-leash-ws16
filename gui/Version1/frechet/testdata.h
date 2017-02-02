#ifndef TESTDATA_H
#define TESTDATA_H

#include <QQuickItem>
#include <QWidget>
#include <QtCore/QObject>
#include <QtCharts/QAbstractSeries>
#include <QObject>


QT_BEGIN_NAMESPACE
class QQuickView;
QT_END_NAMESPACE

QT_CHARTS_USE_NAMESPACE

class Testdata: public QObject
{
    Q_OBJECT
public:
    explicit Testdata(QQuickView *appViewer, QObject *parent = 0);
    //~Testdata();

public slots:
    void cleanDataFile();
    void generateLineTyp(int type, int rowCount, int colCount);
    void update(QAbstractSeries *series);
    void printPointSeries(QAbstractSeries *series);

private:
    QQuickView *m_appViewer;
    QList<QVector<QPointF> > m_data;
    int m_index;
};

#endif // TESTDATA_H


