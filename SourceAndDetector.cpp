#include "SourceAndDetector.h"
#include "ui_SourceAndDetector.h"
#include "Materials.h"

#include <QMessageBox>
#include <QDir>
#include <QFileDialog>
#include <QFile>
// #include <QVector>
#include <QDebug>


SourceAndDetector::SourceAndDetector(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::SourceAndDetector)
{
    ui->setupUi(this);

    ui->lineEdit_XEnergy->setEnabled(false);
    ui->lineEdit_energySpectrum->setEnabled(false);
    ui->pushButton_energySpectrum->setEnabled(false);
    ui->lineEdit_opticalMag->setDisabled(false);

    MaterialsList();
    ui->comboBox_detectorMaterial->addItems(m_materialsList);

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
        QMessageBox::warning(this, "WARNING", "The optical zoom has been set to 1 in the parallel beam mode. Please insert the correct optical zoom.");

    } else {
        return (ui->lineEdit_opticalMag->text().toDouble());
    }
}

void SourceAndDetector::on_radioButton_coneBeam_toggled(bool checked)
{
    m_coneBeam = checked;
}

void SourceAndDetector::on_radioButton_Monochrom_toggled(bool checked)
{
    ui->lineEdit_XEnergy->setEnabled(checked);
        m_monochrome = checked;
}

void SourceAndDetector::on_radioButton_Polychrom_toggled(bool checked)
{
    ui->lineEdit_energySpectrum->setEnabled(checked);
    ui->pushButton_energySpectrum->setEnabled(checked);
    m_polychrome = checked;
}

void SourceAndDetector::getXrayEnergy(std::vector<double>& energyVector, std::vector<double>& spectrumVector)
{
    if (m_monochrome){

        if (!(ui->lineEdit_XEnergy->text().isEmpty())){
            energyVector.push_back((ui->lineEdit_XEnergy->text().toDouble()));
            spectrumVector.push_back(1);
        } else {
            QMessageBox::warning(this, "Error", "Please insert the X-ray energy.");
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
                qDebug() << "energy : " << energy << "\tweight : " << weight;
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
    return (ui->lineEdit_detectorThickness->text().toDouble());
}

double SourceAndDetector::getPixelSizeMM()
{
    if ( (!(ui->lineEdit_pixelSize->text().isEmpty())) && ((ui->lineEdit_pixelSize->text().toDouble()) > 0) ){
        return (ui->lineEdit_pixelSize->text().toDouble());

    } else {

        QMessageBox::warning(this, "ERROR", "Please insert a valid pixel size.");
        return 0;
    }

}

double SourceAndDetector::getFOVmm()
{
    if (!(ui->lineEdit_FOV->text().isEmpty())){
        return (ui->lineEdit_FOV->text().toDouble());

    } else {
        QMessageBox::warning(this, "ERROR", "Please insert a valid field-of-view's size.");
        return 0;
    }
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

// TODO
// std::vector<double> SourceAndDetector::getDetectorEnergyResponse()
// {

// }


