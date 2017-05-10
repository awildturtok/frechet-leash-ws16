import QtQuick 2.5
import QtQuick.Window 2.2
//import QtDataVisualization 1.2
import QtCharts 2.0
import QtQuick 2.0
import QtQuick.Controls 1.0
import QtQuick.Controls.Styles 1.0


ChartView {
    id: chart
    objectName: "chart1"
    title: "Frechét Chart"
    anchors.fill: parent
    antialiasing: true

    //define signals to communicate with cpp classes
    signal sendPoints(string graph, string xCoordinate, string yCoordinate)
    signal curveSignal()
    signal deleteSignal()

    property int curve:1 //1:red 0:blue

    onCurveSignal: {
        switch(curve) {
            case 0: this.curve = 1; break;
            case 1: this.curve = 0; break;
            default:break;
        }
            return "some return value"
    }
    onDeleteSignal: {
        if(curve == 0){//blue
            seriesBLUE.clear();
            scatterBLUE.clear();
        }else if(curve == 1){//red
            seriesRED.clear();
            scatterRED.clear();
        }
            return "some return value"
    }

    LineSeries {
        //graph one
        id: seriesBLUE
        name: "blue"
        color: "steelblue"

        axisX: valueAxisX
        axisY: valueAxisY
    }

    LineSeries {
        //graph two
        id: seriesRED
        name: "red"
        color: "firebrick"

        axisX: valueAxisX
        axisY: valueAxisY
    }

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
        //name: "test"
        color: "#00FF11FF"
        borderColor: "#ff0039A5"
        borderWidth: 0
        axisX: valueAxisX
        axisY: valueAxisY
        upperSeries: LineSeries {
            XYPoint { x: valueAxisX.min; y: valueAxisY.max}
            XYPoint { x: valueAxisX.max; y: valueAxisY.max}
            }
        onClicked: {
            console.log("onClicked: " + point.x + ", " + point.y)
            if(curve == 0){ //state of toggle button //blue
                //send pointinformation over signal to cpp class                
                chart.sendPoints("blue",point.x.toString(),point.y.toString())
                seriesBLUE.append(point.x,point.y)
                //scatterBLUE.append(point.x, point.y)
                scatterBLUE.insert(scatterBLUE.index,point.x,point.y)
                scatterBLUE.index++
                console.log("scatterBLUE point: "+scatterBLUE.at(0));

            }else if (curve == 1){// state of toggle button // red
                chart.sendPoints("red",point.x.toString(),point.y.toString())
                scatterRED.insert(scatterRED.index,point.x,point.y)
                scatterRED.index++
                console.log("scatterRED point: "+scatterRED.at(0));
                seriesRED.append(point.x,point.y);
            }
        }
    }

    ScatterSeries {
        id: scatterRED
        property int index;
        name:"points red"
        axisX: valueAxisX
        axisY: valueAxisY
        color: "lightgreen"
        onHovered: {
            //console.log("onClicked: " + point.x + ", " + point.y+" and test is : "+m.testfunction())
            }
        }

    ScatterSeries {
        property int curve:1
        objectName: "scatterblue"
        id: scatterBLUE
        property int index;
        name:"points blue"
        axisX: valueAxisX
        axisY: valueAxisY
        color: "purple"
        onHovered: {
            //console.log("onClicked: " + point.x + ", " + point.y+" and test is : "+m.testfunction())
        }
    }

}
