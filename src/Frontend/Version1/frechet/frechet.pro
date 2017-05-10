TEMPLATE = app

QT += widgets qml quick charts

CONFIG += c++11

SOURCES += main.cpp \
    mainwindow.cpp

RESOURCES += qml.qrc

HEADERS += \
    mainwindow.h

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

FORMS += \
    mainwindow.ui


