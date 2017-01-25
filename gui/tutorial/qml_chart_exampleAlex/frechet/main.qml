import QtQuick 2.5
import QtQuick.Window 2.2
import QtDataVisualization 1.2
import QtCharts 2.0

Window {
    visible: true
    width: 640
    height: 480
    title: qsTr("Hello World")

    ChartView {
        id: myView
        title: "Line"
        anchors.fill: parent
        antialiasing: true

        LineSeries {
            name: "LineSeries"
            XYPoint { x: 0; y: 0 }
            XYPoint { x: 1.1; y: 2.1 }
            XYPoint { x: 1.9; y: 3.3 }
            XYPoint { x: 2.1; y: 2.1 }
            XYPoint { x: 2.9; y: 4.9 }
            XYPoint { x: 3.4; y: 3.0 }
            XYPoint { x: 4.1; y: 3.3 }
        }
    }

    MouseArea {
        anchors.bottom: myView.bottom
        anchors.top: myView.top
        anchors.right: myView.right
        anchors.left: myView.left
        //anchors.fill: myView
        acceptedButtons: Qt.LeftButton | Qt.RightButton
        onClicked: {
            myView.y
            console.log("x: "+mouseX+", y: "+mouseY)}
    }
}
