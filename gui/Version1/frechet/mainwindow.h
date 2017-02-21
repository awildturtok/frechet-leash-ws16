#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProcess>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void changeSelectedGraph();
    void deleteSelectedGraph();
    void startFrechetCalculation();
    void getPointsFromQML(QString pointInfo);

private:
    Ui::MainWindow *ui;
    QProcess matplotlib;

};

#endif // MAINWINDOW_H
