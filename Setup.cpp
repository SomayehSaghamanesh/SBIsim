#include "Setup.h"
#include "qvalidator.h"
#include "ui_Setup.h"

#include <QMessageBox>
#include <QFileDialog>

Setup::Setup(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Setup)
{
    ui->setupUi(this);

    ui->lineEdit_NumProj->setEnabled(false);

    ui->lineEdit_NumProj->setValidator(new QIntValidator(1, 10000, this));
    ui->lineEdit_SDD->setValidator(new QDoubleValidator(0, 1000000, 6, this));
    ui->lineEdit_SOD->setValidator(new QDoubleValidator(0, 1000000, 6, this));
    ui->lineEdit_SdD->setValidator(new QDoubleValidator(0, 1000000, 6, this));

    connect( ui->lineEdit_NumProj, &QLineEdit::textEdited, this, &Setup::CheckNumProj);
    connect( ui->lineEdit_SDD, &QLineEdit::textEdited, this, &Setup::CheckSDD);
    connect( ui->lineEdit_SOD, &QLineEdit::textEdited, this, &Setup::CheckSOD);
    connect( ui->lineEdit_SdD, &QLineEdit::textEdited, this, &Setup::CheckSdD);

    connect(ui->pushButton_Start, &QPushButton::clicked, this, &Setup::StartButtonClicked);

}

Setup::~Setup()
{
    delete ui;
}

void Setup::on_radioButton_Radiography_toggled(bool checked)
{
    m_radioFlag = checked;
    if (checked){
        ui->lineEdit_NumProj->clear();
    }
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
        return (ui->lineEdit_NumProj->text().toInt());

    } else {

        QMessageBox::warning(this, "ERROR", "Please insert the number of projections");
        return (1);
    }
}

double Setup::getTomoAngleDeg()
{
    return (360/getNumProjections());
}

void Setup::getDistances() // in [m]
{
    if (!(ui->lineEdit_SOD->text().isEmpty()) && ((ui->lineEdit_SOD->text().toDouble()) > 0)){

        m_SOD = 0.001*(ui->lineEdit_SOD->text().toDouble());

    } else {

        QMessageBox::warning(this, "ERROR", "please insert a non-zero source to object distance (mm).");
        return;
    }
    //
    if (!(ui->lineEdit_SdD->text().isEmpty()) && ((ui->lineEdit_SdD->text().toDouble()) > 0.0)){

        m_SdD = 0.001*(ui->lineEdit_SdD->text().toDouble());

    } else {

        QMessageBox::warning(this, "ERROR", "please insert a non-zero source to diffuser distance (mm).");
        return;
    }
    //
    if (!(ui->lineEdit_SDD->text().isEmpty()) && ((ui->lineEdit_SDD->text().toDouble()) > 0.0)){

        m_SDD = 0.001*(ui->lineEdit_SDD->text().toDouble());

    } else {

        QMessageBox::warning(this, "ERROR", "please insert a non-zero source to detector distance (mm).");
        return;
    }

}

void Setup::on_pushButton_ImageDir_clicked()
{
    QString imageDir = QFileDialog::getExistingDirectory(this, tr("Select a directory"), QDir::homePath(), QFileDialog::ShowDirsOnly);
    ui->lineEdit_imageDir->setText(imageDir);
}

QString Setup::getImageDir()
{
    QDir *im_dir = new QDir(ui->lineEdit_imageDir->text());
    if ( (!(ui->lineEdit_imageDir->text().isEmpty()) && (im_dir->exists()))){
        return (ui->lineEdit_imageDir->text());
    } else {
        return (QDir::homePath());
    }

    delete im_dir;
}

void Setup::CheckNumProj()
{
    if (!(ui->lineEdit_NumProj->hasAcceptableInput())){
        QMessageBox::warning(this, "ERROR", "Number of projections must be an integer between 1 and 10000.");
        ui->lineEdit_NumProj->clear();
    }
}

void Setup::CheckSDD()
{
    if (!(ui->lineEdit_SDD->hasAcceptableInput())) {
        QMessageBox::warning(this, "ERROR", "Distance must be a double number between 0 and 1000000.");
        ui->lineEdit_SDD->clear();
    }
}

void Setup::CheckSOD()
{
    if (!(ui->lineEdit_SOD->hasAcceptableInput())) {
        QMessageBox::warning(this, "ERROR", "Distance must be a double number between 0 and 1000000.");
        ui->lineEdit_SOD->clear();
    }
}

void Setup::CheckSdD()
{
    if (!(ui->lineEdit_SdD->hasAcceptableInput())) {
        QMessageBox::warning(this, "ERROR", "Distance must be a double number between 0 and 1000000.");
        ui->lineEdit_SdD->clear();
    }
}
