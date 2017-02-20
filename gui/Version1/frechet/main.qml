import QtQuick 2.5
import QtQuick.Window 2.2
import QtDataVisualization 1.2
import QtCharts 2.0
import QtQuick 2.0
import QtQuick.Controls 1.0
import QtQuick.Controls.Styles 1.0


     ChartView {
        id: chart
        title: "Line"
        anchors.fill: parent
        antialiasing: true
        property int curve:1

        LineSeries {
            //graph two
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
                if(curve == 0){ //state of toggle button //blue
                    seriesBLUE.append(point.x,point.y)
                    //scatterBLUE.append(point.x, point.y)
                    scatterBLUE.insert(scatterBLUE.index,point.x,point.y)
                    scatterBLUE.index++
                    console.log("scatterBLUE point: "+scatterBLUE.at(0));

                }else if (curve == 1){// state of toggle button // red
                    //scatterRED.append(point.x, point.y)
                    console.log("ick hab nen roten punkt")
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
            name:"fancyScatter2"
            axisX: valueAxisX
            axisY: valueAxisY
            color: "firebrick"
            onHovered: {
               //console.log("onClicked: " + point.x + ", " + point.y+" and test is : "+m.testfunction())
            }
        }

        ScatterSeries {
            id: scatterBLUE
            property int index;
            name:"fancyScatter"
            axisX: valueAxisX
            axisY: valueAxisY
            color: "steelblue"
            onHovered: {
               //console.log("onClicked: " + point.x + ", " + point.y+" and test is : "+m.testfunction())
            }
    }
 //Buttons

    //remove button
    Rectangle {
        id: rectangleRemove
        x: 745
        y: 148
        width: 108
        height: 26
        color: "#000000"
    }

    MouseArea {
        id: mouseAreaRemove
        x: 745
        y: 148
        width: 108
        height: 26
//        onClicked: {console.log("onClicked:button " + curve );
//            console.log(rectangleRemove.color);
//            if(rectangleToggleR.color == "#0000ff"){//blue
//                seriesBLUE.clear();
//            }else if(rectangleToggleL.color == "#ff0000"){//red
//                seriesRED.clear();
//            }
//        }
		//test for using cpp methods
        onDoubleClicked: {
            testData.cleanDataFile();
            testData.printPointSeries(chart.series(0));
            testData.printPointSeries(chart.series(1));
		}
    }
    //toggle button
    MouseArea {
        id: mouseAreaToggle
        x: 745
        y: 116
        width: 108
        height: 26
        onClicked: {curve = 2; console.log("onClicked:button " + curve );
            console.log(rectangleToggleR.color);
            rectangleRemove.color = "#000000"; // black
            if(rectangleToggleR.color == "#0000ff"){ //blue
                rectangleToggleR.color = "black";
                rectangleToggleL.color = "red";
                console.log("blue");
                curve = 1;
            }else if(rectangleToggleL.color == "#ff0000"){ //red
                rectangleToggleL.color = "black";
                rectangleToggleR.color = "blue";
                console.log("red");
                curve = 0;
            }
        }
    }

    Rectangle {
        id: rectangleToggleBlack
        x: 773
        y: 116
        width: 50
        height: 26
        color: "black"
    }

    Rectangle {
        id: rectangleToggleR
        x: 823
        y: 116
        width: 30
        height: 26
        color: "black"
    }

    Rectangle {
        id: rectangleToggleL
        x: 745
        y: 116
        width: 28
        height: 26
        color: "red"
    }
}


