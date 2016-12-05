import QtQuick 2.5
import QtQuick.Window 2.2
import QtDataVisualization 1.2

Window {
    visible: true
    width: 640
    height: 480
    title: qsTr("Hello World")

    MainForm {
        anchors.fill: parent
        mouseArea.onClicked: {
            console.log(qsTr('Clicked on background. Text: "' + textEdit.text + '"'))
            // get here the coordinates for frechet elypsies
        }

        Scatter3D {
            id: scatter3D1
            x: 12
            y: 58
            width: 612
            height: 388
            Scatter3DSeries {
                ItemModelScatterDataProxy {
                    itemModel: ListModel {
                        ListElement {
                            x: "1"
                            y: "2"
                            z: "3"
                        }

                        ListElement {
                            x: "2"
                            y: "3"
                            z: "4"
                        }

                        ListElement {
                            x: "3"
                            y: "4"
                            z: "1"
                        }
                    }
                    xPosRole: "x"
                    yPosRole: "y"
                    zPosRole: "z"
                }
            }
        }

    }
}
