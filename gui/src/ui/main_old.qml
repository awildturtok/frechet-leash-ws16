import QtQuick 2.7
import QtQuick.Window 2.2


Rectangle {
    width: 250; height: 175

    Text {
        id: helloText
        anchors.verticalCenter: parent.verticalCenter
        anchors.horizontalCenter: parent.horizontalCenter
        text: "Hello World!!!\n Traditional first app using PyQt5"
        horizontalAlignment: Text.AlignHCenter
    }
}