#include <QApplication>
#include <QQmlApplicationEngine>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QQmlApplicationEngine engine;
    engine.load(QUrl(QStringLiteral("qrc:/main.qml")));

    std::string filename = "test.py";
    std::string command = "python ";
    command += filename;
    system(command.c_str());

    return app.exec();
}

