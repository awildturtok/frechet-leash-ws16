#include <QQuickView>
#include <QtQml/QQmlContext>
#include <QtQml/QQmlEngine>
#include <QtCharts>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <testdata.h>
#include <datahandling.h>
#include <QProcess>
#include <QString>

  QQuickView *viewer;
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    //connect buttons here
    //ToDo: define functions
    connect(ui->btn_change_graph,SIGNAL(clicked(bool)),this,SLOT(changeSelectedGraph()));
    connect(ui->btn_deleteGraph, SIGNAL(clicked(bool)), this, SLOT(deleteSelectedGraph()));
    connect(ui->btn_startCalc, SIGNAL(clicked(bool)), this, SLOT(startFrechetCalculation()));
    viewer = new QQuickView();
    QWidget *container = QWidget::createWindowContainer(viewer,this);
    viewer->setSource(QUrl("qrc:/main.qml"));


    connect(ui->btn_change_graph, SIGNAL(clicked(bool)),viewer->rootObject(), SIGNAL(curveSignal()));
    connect(ui->btn_deleteGraph, SIGNAL(clicked(bool)),viewer->rootObject(), SIGNAL(deleteSignal()));


    //Testdata testData(viewer);
    //ToTest: rootObject returns null pointer
    //connect(viewer->rootObject(), SIGNAL(sendPoints(QString)), this, SLOT(getPointsFromQML(QString)));

    //viewer->rootContext()->setContextProperty("pointData", &pointData);

    ui->verticalLayout_3->addWidget(container);

}

MainWindow::~MainWindow()
{
    delete ui;
}

//change actual graph from blue to red or otherhand
void MainWindow::changeSelectedGraph()
{
}

//delete selected graph
void MainWindow::deleteSelectedGraph(){}


//get information from line edit or from qml chart
void MainWindow::startFrechetCalculation(){

    matplotlib.start("python", QStringList() << ".\\ellipses\\Input.py");
    QString input1= ui->lineEdit_graph_blue->text();
    QString input2= ui->lineEdit_graph_red->text();

    matplotlib.write(input1.toStdString().c_str());
    matplotlib.write("\r\n");
    matplotlib.write(input2.toStdString().c_str());
    matplotlib.write("\r\n");

}

void MainWindow::getPointsFromQML(QString pointInformation){
    qDebug()<< "Points from QML signal" << pointInformation;

    QAbstractSeries *chartView = viewer->rootObject()->findChild<QAbstractSeries*>("scatterblue");
    if(chartView != nullptr){
        QVariant var = chartView->property("curve");
        var.setValue(1);
    }
}
