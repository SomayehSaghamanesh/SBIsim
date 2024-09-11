#include "Setup.h"
#include "ui_Setup.h"

#include <QMessageBox>

Setup::Setup(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Setup)
{
    ui->setupUi(this);

    ui->lineEdit_NumProj->setEnabled(false);

    connect(ui->pushButton_Start, &QPushButton::clicked, this, &Setup::StartButtonClicked);
}

Setup::~Setup()
{
    delete ui;
}

void Setup::on_radioButton_Radiography_toggled(bool checked)
{
    m_radioFlag = checked;
}

void Setup::on_radioButton_Tomography_toggled(bool checked)
{
    ui->lineEdit_NumProj->setEnabled(checked);
    m_tomoFlag = checked;
}

int Setup::getNumProjections()
{
    if (!(ui->lineEdit_NumProj->text().isEmpty()))
    {
        return (ui->lineEdit_NumProj->text().isEmpty());

    } else {

        QMessageBox::warning(this, "ERROR", "Please insert the number of projections");
        return 1;
    }
}

double Setup::getTomoAngleDeg()
{
    return (360/getNumProjections());
}

void Setup::getDistances()
{
    if (!(ui->lineEdit_SOD->text().isEmpty()) && ((ui->lineEdit_SOD->text().toDouble()) > 0.0)){

        m_SOD = (ui->lineEdit_SOD->text().toDouble());

    } else {

        QMessageBox::warning(this, "ERROR", "please insert a valid source to object distance (mm).");
        return;
    }
    //
    if (!(ui->lineEdit_SdD->text().isEmpty()) && ((ui->lineEdit_SdD->text().toDouble()) > 0.0)){

        m_SdD = (ui->lineEdit_SdD->text().toDouble());

    } else {

        QMessageBox::warning(this, "ERROR", "please insert a valid source to diffuser distance (mm).");
        return;
    }
    //
    if (!(ui->lineEdit_SDD->text().isEmpty()) && ((ui->lineEdit_SDD->text().toDouble()) > 0.0)){

        m_SDD = (ui->lineEdit_SDD->text().toDouble());

    } else {

        QMessageBox::warning(this, "ERROR", "please insert a valid source to detector distance (mm).");
        return;
    }

}

// double Setup::getSourceToObjectDist()
// {
//     return (ui->lineEdit_SOD->text().toDouble());
// }

// double Setup::getSourcetoDiffuserDist()
// {
//     return (ui->lineEdit_SdD->text().toDouble());
// }

// double Setup::getSourcetoDetectorDist()
// {
//     return (ui->lineEdit_SDD->text().toDouble());
// }

