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
    void startFrechetCalculation();
    void getPointsFromQML(QString graph, QString xCoord, QString yCoord);
    void changebuttonColor();

private:
    Ui::MainWindow *ui;
    QProcess matplotlib;
    void updateEditLines(QString colour, QString xVal, QString yVal);
    QString color;
};

#endif // MAINWINDOW_H
