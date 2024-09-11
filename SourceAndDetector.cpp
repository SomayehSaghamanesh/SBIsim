#include "SourceAndDetector.h"
#include "ui_SourceAndDetector.h"
#include "Materials.h"

#include <QMessageBox>
#include <QFileDialog>
#include <QFile>
#include <QVector>
#include <QDebug>


SourceAndDetector::SourceAndDetector(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::SourceAndDetector)
{
    ui->setupUi(this);

    ui->lineEdit_XEnergy->setEnabled(false);
    ui->lineEdit_energySpectrum->setEnabled(false);
    ui->pushButton_energySpectrum->setEnabled(false);

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

// double SourceAndDetector::getXrayEnergy()
// {
//     if (!(ui->lineEdit_XEnergy->text().isEmpty())){

//         return (ui->lineEdit_XEnergy->text().toDouble());

//     } else {

//         QMessageBox::warning(this, "ERROR", "Please insert the X-ray energy.");
//         return 0;
//     }
// }

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

        if (!(ui->lineEdit_energySpectrum->text().isEmpty())){

            QString energySpectrum_filename = ui->lineEdit_energySpectrum->text();
            QFile spectrumFile(energySpectrum_filename);
            if (!spectrumFile.open(QFile::ReadOnly | QFile::Text)){// check exist instead of open
                QMessageBox::warning(this, "ERROR", "Energy spectrum pectrum file can NOT open");
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
        } else {
            QMessageBox::warning(this, "ERROR", "Please select a valid energy spectrum file.");
            return;
        }
    } else {
        QMessageBox::warning(this, "ERROR", "Nieghter monochromatic nor polychromatic X-ray has been selected!");
        return;
    }
}

void SourceAndDetector::on_pushButton_energySpectrum_clicked()
{
    QString filter = "Data Files (*.txt, *.dat)";
    QString energySpectrum_filename = QFileDialog::getOpenFileName(this, "Select a file", QDir::homePath(), filter);
    ui->lineEdit_energySpectrum->clear();
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
    return (ui->lineEdit_pixelSize->text().toDouble());
}

double SourceAndDetector::getFOVmm()
{
    return (ui->lineEdit_FOV->text().toDouble());
}

int SourceAndDetector::getNumPixels()
{
    double FOV = getFOVmm();
    int numPixels = FOV/(getPixelSizeMM());

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

double SourceAndDetector::waveLength(double energy)
{
    return (2*pi*energy/hc);
}



