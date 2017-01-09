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

//        LineSeries {
//            //graph one
//            id: series
//            name: "red"
//            color: "red"
//            width: 5
//            XYPoint { x: 0; y: 0 }
//            XYPoint { x: 1.1; y: 2.1 }
//            onClicked: {console.log("onClicked: " + point.x + ", " + point.y);}

//        }

//        LineSeries {
//            //graph two
//            id: series2
//            name: "blue"
//            color: "blue"
//            XYPoint { x: 1; y: 0 }
//            XYPoint { x: 2; y: 0 }
//            XYPoint { x: 3; y: 0 }
//            XYPoint { x: 4; y: 0 }

//            //onClicked: console.log("onClicked: " + point.x + ", " + point.y); //without functionality?!
//        }

        ValueAxis {
               id: valueAxisX
               min: 0
               max: 10
               tickCount: 21
               labelFormat: "%.1f"
           }

        ValueAxis {
               id: valueAxisY
               min: 0
               max: 10
               tickCount: 21
               labelFormat: "%.1f"
           }

        AreaSeries {
            id:backgroundSeries
            name: "test"
            color: "#00FF11FF"
            borderColor: "#ff0039A5"
            borderWidth: 0
            axisX: valueAxisX
            axisY: valueAxisY
            upperSeries: LineSeries {
                XYPoint { x: valueAxisX.min; y: valueAxisY.max}
                XYPoint { x: valueAxisX.max; y: valueAxisY.max }
            }
            onClicked: {
                console.log("onClicked: " + point.x + ", " + point.y)
            }

}

//        MouseArea {
//            anchors.fill: parent
//            acceptedButtons: Qt.LeftButton | Qt.RightButton
//            propagateComposedEvents: true
//            onClicked: {
//                if(mouseX >= 55 && mouseY >= 427){
//                    console.log(mouseX-55)
//                    console.log(mouseY-848)
//                    console.log("added point")
//                    //chart.axisX(series)
//                }
//            }
//            onWheel: {

//            }
//        }

       // focus: true // focus for zooming with keys
//        Keys.onPressed: {
//            /*
//             * zoom
//             * i -> zoom in
//             * o -> zoom out
//             * left -> go left
//             * right -> go right
//             * up -> go up
//             * down -> go down
//              */
//                if (event.key == Qt.Key_I) {
//                    chart.zoomIn()
//                    console.log("axisX: "+backgroundSeries.axisX.min+" and axisY: "+backgroundSeries.axisY)
//                    //backgroundSeries.upperSeries.append(backgroundSeries.axisX,backgroundSeries.axisY)
//                    //backgroundSeries.upperSeries.append()
//                }else if(event.key == Qt.Key_O) {
//                    chart.zoomOut()
//                }else if(event.key == Qt.Key_Left) {
//                    chart.scrollLeft(10)
//                }else if(event.key == Qt.Key_Right) {
//                    chart.scrollRight(10)
//                }else if(event.key == Qt.Key_Up) {
//                    chart.scrollUp(10)
//                }else if(event.key == Qt.Key_Down) {
//                    chart.scrollDown(10)
//                }
//            }

    }
}
