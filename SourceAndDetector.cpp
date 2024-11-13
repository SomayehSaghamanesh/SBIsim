#include "SourceAndDetector.h"
#include "ui_SourceAndDetector.h"
#include "Materials.h"

#include <QMessageBox>
#include <QDir>
#include <QFileDialog>
#include <QFile>
#include <QDebug>


SourceAndDetector::SourceAndDetector(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::SourceAndDetector)
{
    ui->setupUi(this);

    ui->lineEdit_XEnergy->setEnabled(false);
    ui->lineEdit_energySpectrum->setEnabled(false);
    ui->pushButton_energySpectrum->setEnabled(false);
    ui->lineEdit_opticalMag->setEnabled(false);

    MaterialsList();
    ui->comboBox_detectorMaterial->addItems(m_materialsList);

    ui->lineEdit_opticalMag->setValidator(new QDoubleValidator(1, 10000, 6, this));
    ui->lineEdit_XEnergy->setValidator(new QDoubleValidator(0, 10000, 6, this));
    ui->lineEdit_pixelSize->setValidator(new QDoubleValidator(0, 1000, 6, this));
    ui->lineEdit_FOV->setValidator(new QDoubleValidator(0, 10000, 6, this));
    ui->lineEdit_detectorThickness->setValidator(new QDoubleValidator(0, 1000, 6, this));
    // const_cast<QDoubleValidator *>(static_cast<const QDoubleValidator *>(ui->lineEdit_detectorThickness->validator()))->setNotation(QDoubleValidator::StandardNotation);

    connect( ui->lineEdit_opticalMag, &QLineEdit::textEdited, this, &SourceAndDetector::CheckOpticalMag);
    connect( ui->lineEdit_XEnergy, &QLineEdit::textEdited, this, &SourceAndDetector::CheckXEnergy);
    connect( ui->lineEdit_pixelSize, &QLineEdit::textEdited, this, &SourceAndDetector::CheckPixelSize);
    connect( ui->lineEdit_FOV, &QLineEdit::textEdited, this, &SourceAndDetector::CheckOFOV);
    connect( ui->lineEdit_detectorThickness, &QLineEdit::textEdited, this, &SourceAndDetector::CheckDetectorThickness);

}

SourceAndDetector::~SourceAndDetector()
{
    delete ui;
}

void SourceAndDetector::on_radioButton_parBeam_toggled(bool checked)
{
    m_parBeam = checked;
    ui->lineEdit_opticalMag->setEnabled(checked);
}

double SourceAndDetector::getOpticalMag()
{
    if (ui->lineEdit_opticalMag->text().isEmpty()){
        return 1;
        QMessageBox::warning(this, "WARNING", "The optical zoom has been set to 1 in the parallel beam mode. If this is not intended, please insert the correct optical zoom.");

    } else {
        return (ui->lineEdit_opticalMag->text().toDouble());
    }
}

void SourceAndDetector::on_radioButton_coneBeam_toggled(bool checked)
{
    m_coneBeam = checked;
    if (checked){
        ui->lineEdit_opticalMag->clear();
        ui->lineEdit_opticalMag->setEnabled(false);
    }
}

void SourceAndDetector::on_radioButton_Monochrom_toggled(bool checked)
{
    ui->lineEdit_XEnergy->setEnabled(checked);
        m_monochrome = checked;
    if (checked){
        ui->lineEdit_energySpectrum->clear();
        ui->lineEdit_energySpectrum->setEnabled(false);
        ui->pushButton_energySpectrum->setEnabled(false);
    }
}

void SourceAndDetector::on_radioButton_Polychrom_toggled(bool checked)
{
    ui->lineEdit_energySpectrum->setEnabled(checked);
    ui->pushButton_energySpectrum->setEnabled(checked);
    m_polychrome = checked;

    if (checked){
        ui->lineEdit_XEnergy->clear();
        ui->lineEdit_XEnergy->setEnabled(false);
    }
}

void SourceAndDetector::getXrayEnergy(std::vector<double>& energyVector, std::vector<double>& spectrumVector)
{
    if (m_monochrome){

        if ( (!(ui->lineEdit_XEnergy->text().isEmpty())) && ((ui->lineEdit_XEnergy->text().toDouble()) != 0) ){
            energyVector.push_back((ui->lineEdit_XEnergy->text().toDouble()));
            spectrumVector.push_back(1);
        } else {
            QMessageBox::warning(this, "Error", "Please insert a non-zero X-ray energy.");
            return;
        }

    } else if (m_polychrome){

        // QString energySpectrumFilter = "*.txt, *.dat";
        // QString energySpectrum_filename = QFileDialog::getOpenFileName(this, "Open a file", QDir::homePath(), energySpectrumFilter);
        // ui->lineEdit_energySpectrum->clear();
        // ui->lineEdit_energySpectrum->setText(energySpectrum_filename);

        QString energySpectrum_filename = ui->lineEdit_energySpectrum->text();
        QFile spectrumFile(energySpectrum_filename);
        if (!spectrumFile.open(QFile::ReadOnly | QFile::Text)){// check exist instead of open
            QMessageBox::warning(this, "ERROR", "Energy spectrum file can NOT open.");
            return;
        }

        QString line;
        double energy;
        double weight;
        QTextStream spectrumFile_in(&spectrumFile);
        while (spectrumFile_in.readLineInto(&line, 0)){
            QStringList parts = line.split(' ', Qt::SkipEmptyParts);
            if (parts.size() == 2) {
                energy = parts[0].toDouble();
                weight = parts[1].toDouble();
                energyVector.push_back(energy);
                spectrumVector.push_back(weight);
                qDebug() << "Energy : " << energy << "\tWeight : " << weight;
            } else {
                QMessageBox::warning(this, "ERROR", "Invalid file format.");
                return;
            }
        }
        spectrumFile.close();

        // } else {
        //     QMessageBox::warning(this, "ERROR", "Please select a valid energy spectrum file.");
        //     return;
        // }
    } else {
        QMessageBox::warning(this, "ERROR", "Nieghter monochromatic nor polychromatic X-ray has been selected. Aborting ...");
        return;
    }
}

void SourceAndDetector::on_pushButton_energySpectrum_clicked()
{
    ui->lineEdit_energySpectrum->clear();
    QString filter = "Text Files (*.txt);; Data Files (*.dat)";
    QString energySpectrum_filename = QFileDialog::getOpenFileName(this, "Select a file", QDir::homePath(), filter);
    ui->lineEdit_energySpectrum->setText(energySpectrum_filename);
}

QString SourceAndDetector::getDetectorMaterial()
{
    return (ui->comboBox_detectorMaterial->currentText());
}

double SourceAndDetector::getDetectorThickness()
{
    if ( (ui->lineEdit_detectorThickness->text().isEmpty()) || ((ui->lineEdit_detectorThickness->text().toDouble()) == 0) ) {
        QMessageBox::warning(this, "ERROR", "Please insert a non-zero detector thickness.");
    }

    return (ui->lineEdit_detectorThickness->text().toDouble());
}

double SourceAndDetector::getPixelSizeMM()
{
    if ( (ui->lineEdit_pixelSize->text().isEmpty()) || ((ui->lineEdit_pixelSize->text().toDouble()) == 0) ) {
        QMessageBox::warning(this, "ERROR", "Please insert a non-zero detector pixel size.");
    }

    return (ui->lineEdit_pixelSize->text().toDouble());
}

double SourceAndDetector::getFOVmm()
{
    if ( (ui->lineEdit_FOV->text().isEmpty()) || ((ui->lineEdit_FOV->text().toDouble()) == 0) ) {
        QMessageBox::warning(this, "ERROR", "Please insert a valid size for the field of view.");
    }

    return (ui->lineEdit_FOV->text().toDouble());
}

int SourceAndDetector::getNumPixels()
{
    double FOV = getFOVmm();
    double pixelSize = getPixelSizeMM();
    int numPixels = FOV/pixelSize;
    return (numPixels);
}

void SourceAndDetector::MaterialsList()
{
    Materials mat;
    for (int i = 0 ; i < mat.m_properties.size() ; i++)
    {
        m_materialsList.append(mat.m_properties[i].name);
    }
    m_materialsList.sort(Qt::CaseInsensitive);
}

std::vector<double> SourceAndDetector::getWaveNumber(std::vector<double>& energy)
{
    std::vector<double> k; // [1/m]
    for (size_t i = 0 ; i < energy.size() ; i++)
    {
        k.push_back(2*pi*energy[i]/hc);
    }
        return k;
}

void SourceAndDetector::Meshgrid(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y, int& numPixels, double& pixelSize)
{
    // Resize X and Y to be n x n matrices
    X.clear();
    Y.clear();
    X.assign(numPixels, std::vector<double>(numPixels, 0.0));
    Y.assign(numPixels, std::vector<double>(numPixels, 0.0));

    // Fill X and Y coordinate matrices
    for (int i = 0; i < numPixels ; i++) {
        for (int j = 0; j < numPixels; j++) {
            X[i][j] = static_cast<float>((-floor(numPixels/2) + i)*pixelSize*0.001);  // opposite to: X coordinates vary along the columns
            Y[i][j] = static_cast<float>((-floor(numPixels/2) + j)*pixelSize*0.001);  // opposite to: Y coordinates vary along the rows
        }
    }
}

void SourceAndDetector::DetectorCoordinates(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y, std::vector<std::vector<double>>& rsqr, int& numPixels)
{
    rsqr.clear();
    rsqr.assign(numPixels, std::vector<double>(numPixels, 0.0));

    for (int i = 0 ; i < numPixels ; i++){
        for (int j = 0 ; j < numPixels ; j++){
            // r*r in m
            rsqr[i][j] = std::pow(X[i][j], 2) + std::pow(Y[i][j], 2); // [m^2]
        }
    }
}

void SourceAndDetector::CheckOpticalMag()
{
    if (!(ui->lineEdit_opticalMag->hasAcceptableInput())){
        QMessageBox::warning(this, "ERROR", "Optical magnification factor must be a double number greater than 1.");
        ui->lineEdit_opticalMag->clear();
    }
}

// void SourceAndDetector::CheckValues(QLineEdit *lineEdit)
// {
//     if ( (!(lineEdit->hasAcceptableInput())) || ((lineEdit->text().toDouble())==0)){
//         QMessageBox::warning(this, "ERROR", "Please insert a double number between 0 and 10000.");
//         lineEdit->clear();
//     }
// }

void SourceAndDetector::CheckXEnergy()
{
    if (!(ui->lineEdit_XEnergy->hasAcceptableInput())) {
        QMessageBox::warning(this, "ERROR", "X-ray energy must be a double number between 0 and 10000.");
        ui->lineEdit_XEnergy->clear();
    }
}

void SourceAndDetector::CheckPixelSize()
{
    if (!(ui->lineEdit_pixelSize->hasAcceptableInput())) {
        QMessageBox::warning(this, "ERROR", "Pixel size must be a double number between 0 and 1000.");
        ui->lineEdit_pixelSize->clear();
    }
}

void SourceAndDetector::CheckOFOV()
{
    if (!(ui->lineEdit_FOV->hasAcceptableInput())) {
        QMessageBox::warning(this, "ERROR", "Field of view must be a double number between 0 and 10000.");
        ui->lineEdit_FOV->clear();
    }
}

void SourceAndDetector::CheckDetectorThickness()
{
    if (!(ui->lineEdit_detectorThickness->hasAcceptableInput())) {
        QMessageBox::warning(this, "ERROR", "Detector thickness must be a double number between 0 and 1000.");
        ui->lineEdit_detectorThickness->clear();
    }
}

// TODO
// std::vector<double> SourceAndDetector::getDetectorEnergyResponse()
// {

// }


