TEMPLATE = app

QT += widgets qml quick charts

CONFIG += c++11 noKeywords

SOURCES += main.cpp

RESOURCES += qml.qrc
# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

win32:CONFIG(release): LIBS += -L../../../../../../AppData/Local/Programs/Python/Python35-32/libs/ -lpython35
#else:win32:CONFIG(debug, debug|release): LIBS += -L../../../../../../AppData/Local/Programs/Python/Python35-32/libs/ -lpython35d
else:unix: LIBS += -L../../../../../../AppData/Local/Programs/Python/Python35-32/libs/ -lpython35

INCLUDEPATH += ../../../../../../AppData/Local/Programs/Python/Python35-32/libs
DEPENDPATH += ../../../../../../AppData/Local/Programs/Python/Python35-32/libs
