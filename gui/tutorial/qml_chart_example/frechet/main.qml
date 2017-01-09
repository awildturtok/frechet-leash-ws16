import QtQuick 2.5
import QtQuick.Window 2.2
import QtDataVisualization 1.2
import QtCharts 2.0

Window {
    visible: true
    width: 900
    height: 900
    title: qsTr("Hello World")

    ChartView {
        id: chart
        title: "Line"
        anchors.fill: parent
        antialiasing: true

        LineSeries {
            //graph one
            id: series
            name: "red"
            color: "red"
            XYPoint { x: 1; y: 1 }
            XYPoint { x: 2; y: 2 }
            XYPoint { x: 3; y: 3 }
            XYPoint { x: 4; y: 4 }

            //onClicked: console.log("onClicked: " + point.x + ", " + point.y); //without functionality?!
        }

        LineSeries {
            //graph two
            id: series2
            name: "blue"
            color: "blue"
            XYPoint { x: 1; y: 0 }
            XYPoint { x: 2; y: 0 }
            XYPoint { x: 3; y: 0 }
            XYPoint { x: 4; y: 0 }

            //onClicked: console.log("onClicked: " + point.x + ", " + point.y); //without functionality?!
        }

        MouseArea {
            anchors.fill: parent
            acceptedButtons: Qt.LeftButton | Qt.RightButton
            onClicked: {
                if(mouseX >= 55 && mouseY >= 427){
                    console.log(mouseX-55)
                    console.log(mouseY-848)
                    series.append((mouseX-55),(mouseY-427))
                    chart.axisX(series)
                }
            }
        }

        focus: true // focus for zooming with keys
        Keys.onPressed: {
            /*
             * zoom
             * i -> zoom in
             * o -> zoom out
             * left -> go left
             * right -> go right
             * up -> go up
             * down -> go down
              */
                if (event.key == Qt.Key_I) {
                    chart.zoom(1.3)
                }else if(event.key == Qt.Key_O) {
                    chart.zoom(0.1)
                }else if(event.key == Qt.Key_Left) {
                    chart.scrollLeft(10)
                }else if(event.key == Qt.Key_Right) {
                    chart.scrollRight(10)
                }else if(event.key == Qt.Key_Up) {
                    chart.scrollUp(10)
                }else if(event.key == Qt.Key_Down) {
                    chart.scrollDown(10)
                }
            }
    }
}
