#include "DiffuserAndObject.h"
#include "ui_DiffuserAndObject.h"
#include "Materials.h"

#include <QMessageBox>
#include <QDebug>
#include <QFile>
#include <QFileDialog>
#include <algorithm>


DiffuserAndObject::DiffuserAndObject(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::DiffuserAndObject)
{
    ui->setupUi(this);

    ui->lineEdit_CDiameter->setEnabled(false);
    ui->lineEdit_CHeight->setEnabled(false);
    ui->lineEdit_SDiameter->setEnabled(false);

    MaterialsList();
    ui->comboBox_OMaterial->addItems(m_materialsList);
    ui->comboBox_numMVSlices->addItems(m_num_MVS_slices);
    ui->comboBox_GritMaterial->addItems(m_materialsList);
    ui->comboBox_BaseMaterial->addItems(m_materialsList);

    // MagFactors();

}

DiffuserAndObject::~DiffuserAndObject()
{
    delete ui;
}

void DiffuserAndObject::on_radioButton_Cylinder_toggled(bool checked)
{
    ui->lineEdit_CDiameter->setEnabled(checked);
    ui->lineEdit_CHeight->setEnabled(checked);
    m_cylinder = checked;
}

void DiffuserAndObject::on_radioButton_Sphere_toggled(bool checked)
{
    ui->lineEdit_SDiameter->setEnabled(checked);
    m_sphere = checked;
}

void DiffuserAndObject::on_radioButton_VObject_toggled(bool checked)
{
    ui->lineEdit_VObject->setEnabled(checked);
    ui->pushButton_VObject->setEnabled(checked);
    m_virtualObject = checked;
}

void DiffuserAndObject::on_pushButton_VObject_clicked()
{
    QString filter = "Image Files (*.raw, *.img)";
    QString vObject_filename = QFileDialog::getOpenFileName(this, "Select a file", QDir::homePath(), filter);
    ui->lineEdit_VObject->clear();
    ui->lineEdit_VObject->setText(vObject_filename);
}

double DiffuserAndObject::getCylinderDiameter()
{
    return (ui->lineEdit_CDiameter->text().toDouble());
}

double DiffuserAndObject::getCylinderHeight()
{
    return (ui->lineEdit_CHeight->text().toDouble());
}

double DiffuserAndObject::getSphereDiameter()
{
    return (ui->lineEdit_SDiameter->text().toDouble());
}

QVector<QVector<QVector<float>>> DiffuserAndObject::getVirtualObject()
{
    if (m_virtualObject){

        QString vObject_filename = ui->lineEdit_VObject->text();
        QFile VObjectFile(vObject_filename);
        if (!VObjectFile.open(QFile::ReadOnly | QFile::Text)){// check exist instead of open
            QMessageBox::warning(this, "ERROR", "Image file can NOT open");
            return {};
        }
            // TODO
        // QString line;
        // float energy;
        // float weight;
        // QTextStream spectrumFile_in(&spectrumFile);
        // while (spectrumFile_in.readLineInto(&line, 0)){
        //     QStringList parts = line.split(' ', Qt::SkipEmptyParts);
        //     if (parts.size() == 2) {
        //         energy = parts[0].toFloat();
        //         weight = parts[1].toFloat();
        //         m_energyVector.push_back(energy);
        //         m_SpectrumVector.push_back(weight);
        //         qDebug() << "energy : " << energy << "\tweight : " << weight;
        //     } else {
        //         QMessageBox::warning(this, "ERROR", "Invalid file format.");
        //         return;
        //     }
        VObjectFile.close();
        return m_vObject;

    } else {
        QMessageBox::warning(this, "ERROR", "Virtual Object NOT found!");
        return {};
    }
}

void DiffuserAndObject::MaterialsList()
{
    Materials mat;
    for (int i = 0 ; i < mat.m_properties.size() ; i++)
    {
        m_materialsList.append(mat.m_properties[i].name);
    }
    m_materialsList.sort(Qt::CaseInsensitive);
}

QString DiffuserAndObject::getObjectMaterial()
{
    return (ui->comboBox_OMaterial->currentText());
}

int DiffuserAndObject::getNumMVSlices(double thickness, int& numVoxelsInZ, double pixelSize)
{
    numVoxelsInZ = static_cast<int>(thickness/(pixelSize));
    int numSlices = 1; // no slice
    // qDebug() << "numVoxelsInZ" << numVoxelsInZ;
    // qDebug() << "voxel_size =" << pixelSize;

    if ((ui->comboBox_numMVSlices->currentText()) == "large"){
        QMessageBox::warning(this, "WARNING", "Single-voxel slices have been selected.");
        numSlices = numVoxelsInZ;// 1 slice = 1 voxel

    } else if ((ui->comboBox_numMVSlices->currentText()) == "medium"){
        numSlices = std::max((numVoxelsInZ/2), 1); // 1 slice = 2 voxels

    } else if ((ui->comboBox_numMVSlices->currentText()) == "few"){
        numSlices = std::max((numVoxelsInZ/100*20), 1); // 20% of total voxels along thickness

    } else {
        QMessageBox::warning(this, "WARNING", "only one slice has been selected for MVS process.");
        // return 1; // no slice
    }

        return numSlices;
}

int DiffuserAndObject::getNumInterpolations()
{
    if ( (!(ui->lineEdit_NumInterp->text().isEmpty())) && ((ui->lineEdit_NumInterp->text().toInt()) > 1) ) {
        return (ui->lineEdit_NumInterp->text().toInt());
    } else {
        QMessageBox::warning(this, "WARNING", "One interpolation was selected by user.");
        return 1;
    }
}

double DiffuserAndObject::getGritSizeMM()
{
    return (ui->lineEdit_gritSize->text().toDouble());
}

QString DiffuserAndObject::getGritMaterial()
{
    return (ui->comboBox_GritMaterial->currentText());
}

QString DiffuserAndObject::getBaseMaterial()
{
    return (ui->comboBox_BaseMaterial->currentText());
}

double DiffuserAndObject::getObjectThickness()
{
    if (m_sphere){
        return (getSphereDiameter());

    } else if (m_cylinder){
        return (getCylinderDiameter());
    }
    // } else if (m_vObject){
    //     float thickness = (m_vObject[0][0].size())*getPixelSizeMM()*SOD/SDD; // only if we assume the voxel size in v_object is the same as simulation here
    //         return (thickness);
    // }
    return 0;
}

double DiffuserAndObject::getDiffuserThickness()
{
    if (!(ui->lineEdit_diffThickness->text().isEmpty()) && ((ui->lineEdit_diffThickness->text().toFloat()) > 0)){
        return (ui->lineEdit_diffThickness->text().toDouble());
    } else {
        QMessageBox::warning(this, "ERROR", "Please insert a valid diffuser thickness (mm).");
        return 0;
    }

}

// void DiffuserAndObject::MagFactors()
// {
//     // magnification factors at each slice inside obj/diff
//     m_dM_obj.clear();
//     m_dM_diff.clear();
//     // mag factors based on the fresnel scaling theory
//     m_fM_obj.clear();
//     m_fM_diff.clear();

//     m_objThickness = getObjectThickness();
//     m_diffThickness = getDiffuserThickness();

//     Setup setup;
//     setup.getDistances();

//     SourceAndDetector det;
//     m_pixelSize = det.m_pixelSizeMM;
//     m_numPixels = det.getNumPixels();

//     // global magnification factors
//     m_M_obj = (setup.m_SDD)/(setup.m_SOD);
//     m_M_diff = (setup.m_SDD)/(setup.m_SdD);

//     m_pixelObj = m_pixelSize/m_M_obj;
//     m_pixelDiff = m_pixelSize/m_M_diff;

//     m_objThicknessAsIndex = m_objThickness/m_pixelObj;
//     m_diffThicknessAsIndex = m_diffThickness/m_pixelDiff;

//     for (int i = 0 ; i < m_objThicknessAsIndex ; i++)
//     {
//         m_dM_obj.push_back((setup.m_SDD)/((setup.m_SOD) - m_objThickness/2 + (i - 1)*m_pixelSize/m_M_obj)); // dynamic magnification factors of slices
//         m_fM_obj.push_back(((setup.m_SOD) - m_objThickness/2 + (i - 1)*m_pixelSize/m_M_obj)/((setup.m_SOD) - m_objThickness/2));; // Magnification after applying Fresnel scaling theory
//     }

//     for (int i = 0 ; i < m_diffThicknessAsIndex ; i++)
//     {
//         m_dM_diff.push_back((setup.m_SDD)/((setup.m_SdD) - m_diffThickness/2 + (i - 1)*m_pixelSize/m_M_diff)); // dynamic magnification factors of slices
//         m_fM_diff.push_back(((setup.m_SdD) - m_diffThickness/2 + (i - 1)*m_pixelSize/m_M_diff)/((setup.m_SdD) - m_diffThickness/2)); // Magnification after applying Fresnel scaling theory
//     }
// }

