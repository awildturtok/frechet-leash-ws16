#include <QQuickView>
#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    //connect buttons here
    QQuickView *view = new QQuickView();
    QWidget *container = QWidget::createWindowContainer(view,this);
    view->setSource(QUrl("qrc:/main.qml"));
    ui->gridLayout->addWidget(container);
}

MainWindow::~MainWindow()
{
    delete ui;
}
