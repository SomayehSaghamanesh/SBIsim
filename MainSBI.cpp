#include "MainSBI.h"
#include "ui_MainSBI.h"
#include "Simulation.h"

#include <QApplication>

MainSBI::MainSBI(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainSBI)
{
    ui->setupUi(this);

}

MainSBI::~MainSBI()
{
    delete ui;
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

