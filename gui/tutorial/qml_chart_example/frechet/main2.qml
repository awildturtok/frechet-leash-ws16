
import QtQuick 2.5
import QtQuick.Window 2.2

Window {
    visible: true
    width: 900
    height: 900
    title: qsTr("Hello World")
    property int xpos
    property int ypos
Rectangle {
    width: 360
    height: 360


    Canvas {
        id: myCanvas
        anchors.fill: parent

        onPaint: {
            var ctx = getContext('2d')
            ctx.fillStyle = "red"
            ctx.fillRect(xpos-1, ypos-1, 3, 3)

        }

        MouseArea{
            anchors.fill: parent
            onPressed: {
                xpos = mouseX
                ypos = mouseY
                myCanvas.requestPaint()
                console.log("onClicked: " + point.x + ", " + point.y);
            }
            onMouseXChanged: {
                xpos = mouseX
                ypos = mouseY
                myCanvas.requestPaint()
            }
            onMouseYChanged: {
                xpos = mouseX
                ypos = mouseY
                myCanvas.requestPaint()
            }
        }

    }

}
}
