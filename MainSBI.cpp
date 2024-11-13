#include "MainSBI.h"
#include "ui_MainSBI.h"
#include "Simulation.h"

#include <QApplication>
#include <QPixmap>


MainSBI::MainSBI(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainSBI)
{
    ui->setupUi(this);

    // QPixmap pix("../../img/SBIsim_cover.png");
    // QString style = QString("MainSBI { image: url(:/img/SBIsim_cover.png); image-position: center; }");
    // style.repeated(1);
    // this->setStyleSheet(style);

    // QString relativePath = "../../img/SBIsim_cover.png";
    // QString absolutePath = QDir::current().absoluteFilePath(relativePath);
    // QPixmap pic("/home/somayeh/QT_GUI/SBIsim/img/SBIsim_cover.png");
    // int w = ui->label_cover->width();
    // int h = ui->label_cover->height();
    // ui->label_cover->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
    // ui->label_cover->setDisabled(true);
    // this->setStyleSheet("MainSBI { background-image: url(:/Resources/img/SBIsim_cover.png);"
    //                     "background-size: contain;"
    //                     "background-repeat: no-repeat; background-position: center; }");

}

MainSBI::~MainSBI()
{
    delete ui;
}

// void MainSBI::showEvent(QShowEvent *event) {

//     QMainWindow::showEvent(event);  // Call the base class implementation
//     QString relativePath = "../../img/SBIsim_cover.png";
//     QString absolutePath = QDir::current().absoluteFilePath(relativePath);
//     QPixmap pix(absolutePath);
//     int w = ui->label_cover->width();
//     int h = ui->label_cover->height();
//     ui->label_cover->setPixmap(pix.scaled(w, h, Qt::KeepAspectRatio));
//     ui->label_cover->setDisabled(true);
// }

void MainSBI::paintEvent(QPaintEvent *event) {

    QMainWindow::paintEvent(event);  // Call the base class implementation

    QPixmap pixmap(":/Resources/img/SBIsim_cover.png");
    if (pixmap.isNull()) {
        return; // Make sure the image is loaded
    }

    QSize windowSize = this->size();

    // Calculate the scaling factor to fit the image within the window
    double scaleFactor = std::min(
        static_cast<double>(windowSize.width()) / pixmap.width(),
        static_cast<double>(windowSize.height()) / pixmap.height()
        );

    // Calculate the scaled size based on the scaling factor
    QSize scaledSize = pixmap.size() * scaleFactor * 0.8;

    // Scale the pixmap to fit the window while keeping the aspect ratio
    QPixmap scaledPixmap = pixmap.scaled(scaledSize, Qt::KeepAspectRatio, Qt::SmoothTransformation);

    // Calculate the position to center the image in the window
    int x = (windowSize.width() - scaledPixmap.width()) / 2;
    int y = (windowSize.height() - scaledPixmap.height()) / 2;

    // Draw the scaled image centered in the window
    QPainter painter(this);
    painter.setOpacity(0.6);
    painter.drawPixmap(x, y, scaledPixmap);
}

void MainSBI::on_actionSimulation_triggered()
{
    Simulation *simWindow = new Simulation(this);
    setCentralWidget(simWindow);
    simWindow->show();
}

void MainSBI::on_actionExit_triggered()
{
    QApplication::quit();
}

