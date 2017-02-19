TEMPLATE = app

QT += widgets qml quick charts

CONFIG += c++11

SOURCES += main.cpp \
    testdata.cpp \
    mainwindow.cpp \
    datahandling.cpp

RESOURCES += qml.qrc

HEADERS += \
    testdata.h \
    mainwindow.h \
    datahandling.h

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

FORMS += \
    mainwindow.ui


