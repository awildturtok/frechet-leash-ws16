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
    QQuickView *viewer = new QQuickView();
    QWidget *container = QWidget::createWindowContainer(viewer,this);

    connect(viewer->engine(), &QQmlEngine::quit, viewer, &QApplication::quit);

    Testdata testData(viewer);
    viewer->rootContext()->setContextProperty("testData", &testData);
    viewer->setSource(QUrl("qrc:/main.qml"));
    ui->verticalLayout->addWidget(container);

}

MainWindow::~MainWindow()
{
    delete ui;
}
