#include <QQuickView>
#include <QtQml/QQmlContext>
#include <QtQml/QQmlEngine>
#include <QtCharts>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QProcess>
#include <QString>
#include <sstream>
#include<QWidget>
#include <QQuickItem>

  QQuickView *viewer;
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    connect(ui->btn_startCalc, SIGNAL(clicked(bool)), this, SLOT(startFrechetCalculation()));
    viewer = new QQuickView();
    QWidget *container = QWidget::createWindowContainer(viewer,this);
    viewer->setSource(QUrl("qrc:/main.qml"));

    connect(viewer->rootObject(), SIGNAL(sendPoints(QString, QString, QString)), this, SLOT(getPointsFromQML(QString, QString, QString)));
    connect(ui->btn_change_graph, SIGNAL(clicked(bool)),viewer->rootObject(), SIGNAL(curveSignal()));
    connect(ui->btn_deleteGraph, SIGNAL(clicked(bool)), ui->lineEdit_graph_blue, SLOT(clear()));
    connect(ui->btn_deleteGraph, SIGNAL(clicked(bool)), ui->lineEdit_graph_red, SLOT(clear()));
    connect(ui->btn_deleteGraph, SIGNAL(clicked(bool)),viewer->rootObject(), SIGNAL(deleteSignal()));
    connect(ui->slider_ellipsies,SIGNAL(valueChanged(int)), ui->lb_ellipsiesValue, SLOT(setNum(int)));
    connect(ui->slider_samples, SIGNAL(valueChanged(int)), ui->lb_sampleValue, SLOT(setNum(int)));
    connect(ui->btn_change_graph,SIGNAL(clicked(bool)), this, SLOT(changebuttonColor()));
    ui->btn_change_graph->setStyleSheet("background-color: red");
    this->color = "red";

    ui->verticalLayout_3->addWidget(container);

}

MainWindow::~MainWindow()
{
    delete ui;
}

//get information from line edit or from qml chart
void MainWindow::startFrechetCalculation(){

    //n number of ellypsie contour (height lines) optimal 7, range 3-20
    //l list of special ellypsie contours
    //s number of samples -> resolution , optimal 100, range 30-200
    int numSamples = ui->slider_samples->value();
    int numEllipsies = ui->slider_ellipsies->value();

    matplotlib.start("python", QStringList() << ".\\ellipses\\Input.py" << "-n" << QString::number(numEllipsies) << "-s" << QString::number(numSamples));
    QString input1= ui->lineEdit_graph_blue->text();
    QString input2= ui->lineEdit_graph_red->text();

    matplotlib.write(input1.toStdString().c_str());
    matplotlib.write("\r\n");
    matplotlib.write(input2.toStdString().c_str());
    matplotlib.write("\r\n");

}

void MainWindow::getPointsFromQML(QString graphCol, QString xCoord, QString yCoord){

    updateEditLines(graphCol, xCoord, yCoord);

    QAbstractSeries *chartView = viewer->rootObject()->findChild<QAbstractSeries*>("scatterblue");
    if(chartView != nullptr){
        QVariant var = chartView->property("curve");
        var.setValue(1);

    }
}

void MainWindow::updateEditLines(QString color, QString x, QString y){

    std::stringstream pointsOnRedLine;
    std::stringstream pointsOnBlueLine;
    if (color=="red")
    {
        if(ui->lineEdit_graph_red->text()==""){
            pointsOnRedLine << x.toStdString() << ";" << y.toStdString();
        }
        else{
            QString actualText = ui->lineEdit_graph_red->text();
            pointsOnRedLine << actualText.toStdString() << "," << x.toStdString() << ";" << y.toStdString();
        }
        QString lineRed = QString::fromStdString(pointsOnRedLine.str());
        ui->lineEdit_graph_red->setText(lineRed);
    }else
    {
        if(ui->lineEdit_graph_blue->text()==""){
            pointsOnBlueLine << x.toStdString() << ";" << y.toStdString();
        }
        else{
            QString actualText = ui->lineEdit_graph_blue->text();
            pointsOnBlueLine << actualText.toStdString() << "," << x.toStdString() << ";" << y.toStdString();
        }
        QString lineBlue = QString::fromStdString(pointsOnBlueLine.str());
        ui->lineEdit_graph_blue->setText(lineBlue);
    }
}

void MainWindow::changebuttonColor(){
    if (this->color == "red"){

        ui->btn_change_graph->setStyleSheet("background-color: blue");
        this->color = "blue";
    }
    else
    {
        ui->btn_change_graph->setStyleSheet("background-color: red");
        this->color= "red";
    }

}
