#include <QQuickView>
#include <QtQml/QQmlContext>
#include <QtQml/QQmlEngine>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <testdata.h>


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

    QQuickView *viewer = new QQuickView();
    QWidget *container = QWidget::createWindowContainer(viewer,this);

    connect(viewer->engine(), &QQmlEngine::quit, viewer, &QApplication::quit);

    Testdata testData(viewer);
    viewer->rootContext()->setContextProperty("testData", &testData);
    viewer->setSource(QUrl("qrc:/main.qml"));
    ui->verticalLayout_3->addWidget(container);

}

MainWindow::~MainWindow()
{
    delete ui;
}

//change actual graph from blue to red or otherhand
void MainWindow::changeSelectedGraph(){}

//delete selected graph
void MainWindow::deleteSelectedGraph(){}


//get information from line edit or from qml chart
void MainWindow::startFrechetCalculation(){}
