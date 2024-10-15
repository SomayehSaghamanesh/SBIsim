#include "DiffuserAndObject.h"
#include "ui_DiffuserAndObject.h"
#include "Materials.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"

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
    ui->comboBox_objectMaterial->addItems(m_materialsList);
    ui->comboBox_numMVSlices->addItems(m_num_MVS_slices);
    ui->comboBox_GritMaterial->addItems(m_materialsList);
    ui->comboBox_BaseMaterial->addItems(m_materialsList);
    ui->comboBox_gritDensity->addItems(m_gritDensity);

    ui->lineEdit_NumInterp->setText("1");

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

QString DiffuserAndObject::getObjectName()
{
    if (m_sphere){
        return ("Sphere");

    } else if (m_cylinder){
        return ("Cylinder");

    } else if (m_virtualObject){
        return ("VirtualObject");

    } else {
        QMessageBox::warning(this, "ERROR", "No object was selected. Aborting ...");
        return "";
    }
}

void DiffuserAndObject::on_pushButton_VObject_clicked()
{
    ui->lineEdit_VObject->clear();
    // QString filter = "Image Files (*.raw, *.img)";
    QString vObject_filename = QFileDialog::getOpenFileName(this, "Select a file", QDir::homePath());
    ui->lineEdit_VObject->setText(vObject_filename);
}

double DiffuserAndObject::getCylinderDiameter()
{
    return (0.001*(ui->lineEdit_CDiameter->text().toDouble())); // [m]
}

double DiffuserAndObject::getCylinderHeight()
{
    return (0.001*(ui->lineEdit_CHeight->text().toDouble())); // [m]
}

double DiffuserAndObject::getSphereDiameter()
{
    return (0.001*(ui->lineEdit_SDiameter->text().toDouble())); // [m]
}

double DiffuserAndObject::getVObjectMag()
{
    if ( (m_virtualObject) && !(ui->lineEdit_vObjectMag->text().isEmpty()) ){
       return (ui->lineEdit_vObjectMag->text().toDouble());

    } else {
        QMessageBox::warning(this, "WARNING", "Please insert the virtual object magnification or it will be set to 1.");
        return 1;
    }
}

std::vector<std::vector<std::vector<float>>> DiffuserAndObject::getVirtualObject(const double& mag_vObject,
                                                                const int& numPixels, const int& numPixelsZ, const double& PixelSizeObj)
{
    QString vObject_filename = ui->lineEdit_VObject->text();

    if ( (m_virtualObject) && (QFile::exists(vObject_filename)) && (!(ui->lineEdit_VObject->text().isEmpty())) ){
    } else {
        QMessageBox::warning(this, "ERROR", "Virtual image file can NOT open");
        return {};
    }

    std::vector<std::vector<std::vector<float>>> vObject(numPixels, std::vector<std::vector<float>>(numPixels, std::vector<float>(numPixelsZ)));

    // image type
    using ImageType = itk::Image<float, 3>;
    using ReaderType =itk::ImageFileReader<ImageType>;
    ReaderType::Pointer reader = ReaderType::New();

    // set the file name
    reader->SetFileName(vObject_filename.toStdString());
    try{
        reader->Update();
    } catch (itk::ExceptionObject & err){
        qDebug() << "Error reading the virtual image : " << err.what();
        throw;
    }
    // get the image
    ImageType::Pointer vImage = reader->GetOutput();

    // resampler
    using ResamplerFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
    ResamplerFilterType::Pointer resampleFilter = ResamplerFilterType::New();
    // set the transform
    using TransformType = itk::AffineTransform<double, 3>;
    TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    resampleFilter->SetTransform(transform);
    // set (linear) interpolator
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    resampleFilter->SetInterpolator(interpolator);
    // set the output spacing
    ImageType::SpacingType inputSpacing = (vImage->GetSpacing())/mag_vObject; // the old spacing on the object after demagnification
    ImageType::SpacingType outputSpacing; // real pixel size on the object phantom
    outputSpacing[0] = PixelSizeObj;
    outputSpacing[1] = PixelSizeObj;
    outputSpacing[2] = PixelSizeObj;
    resampleFilter->SetOutputSpacing(outputSpacing);
    // set the size of output image
    ImageType::RegionType region = vImage->GetLargestPossibleRegion();
    ImageType::SizeType inputSize = region.GetSize();
    ImageType::SizeType outputSize;
    outputSize[0] = static_cast<unsigned int>(inputSize[0] * (inputSpacing[0]/outputSpacing[0]));
    outputSize[1] = static_cast<unsigned int>(inputSize[1] * (inputSpacing[1]/outputSpacing[1]));
    outputSize[2] = static_cast<unsigned int>(inputSize[2] * (inputSpacing[2]/outputSpacing[2]));
    resampleFilter->SetSize(outputSize);
    // set the origin and direction as the same as the input image
    resampleFilter->SetOutputOrigin(vImage->GetOrigin());
    resampleFilter->SetOutputDirection(vImage->GetDirection());
    // set input image to the resampler
    resampleFilter->SetInput(vImage);
    // generate the resampled image
    resampleFilter->Update();
    // get final resampled image
    ImageType::Pointer resampledImage = resampleFilter->GetOutput();
    // check if axes order is correct
    bool reverseAxes = false;

    // calculate cropping offsets
    ImageType::RegionType finalRegion = resampledImage->GetLargestPossibleRegion();
    ImageType::SizeType finalSize = finalRegion.GetSize();
    int offsetX = (finalSize[0] - numPixels)/2;
    int offsetY = (finalSize[1] - numPixels)/2;
    int offsetZ = (finalSize[2] - numPixels)/2;
    //
    offsetX = std::max(0, offsetX);
    offsetY = std::max(0, offsetY);
    offsetZ = std::max(0, offsetZ);

    // iterate through the final image
    itk::ImageRegionIterator<ImageType> imageIter(resampledImage, finalRegion);
    int xIndex = 0, yIndex = 0, zIndex = 0;

    // Iterate over the cropped region and assign values to m_vObject
    for (unsigned int z = offsetZ; z < std::min(static_cast<unsigned int>(finalSize[2]), static_cast<unsigned int>(offsetZ + numPixelsZ)); ++z) {
        zIndex = z - offsetZ;
        for (unsigned int y = offsetY; y < std::min(static_cast<unsigned int>(finalSize[1]), static_cast<unsigned int>(offsetY + numPixels)); ++y) {
            yIndex = y - offsetY;
            for (unsigned int x = offsetX; x < std::min(static_cast<unsigned int>(finalSize[0]), static_cast<unsigned int>(offsetX + numPixels)); ++x) {
                xIndex = x - offsetX;

                ImageType::IndexType index = {{x, y, z}};
                float pixelValue = resampledImage->GetPixel(index);

                if (reverseAxes){
                    // store ZYX instead of XYZ
                    vObject[zIndex][yIndex][xIndex] = pixelValue;

                } else {
                    vObject[xIndex][yIndex][zIndex] = pixelValue;
                }
            }
        }
    }

    // while (!imageIter.IsAtEnd()){
    //     ImageType::IndexType index = imageIter.GetIndex();
    //     float pixelValue = imageIter.Get();
    //     m_vObject[index[0]][index[1]][index[2]] = pixelValue;
    //     ++imageIter;
    // }

    return vObject;
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
    return (ui->comboBox_objectMaterial->currentText());
}

int DiffuserAndObject::getNumMVSlices(const double& thickness, int& numVoxelsInZ, const double& pixelSize, const double& magnification)
{
    numVoxelsInZ = static_cast<int>(thickness/(pixelSize/magnification));
    int numSlices = 1; // no slice
    // qDebug() << "numVoxelsInZ" << numVoxelsInZ;
    // qDebug() << "voxel_size =" << pixelSize;

    if ((ui->comboBox_numMVSlices->currentText()) == "large"){
        // QMessageBox::warning(this, "WARNING", "Single-voxel slices have been selected.");
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

double DiffuserAndObject::getGritSize()
{
    return (0.001*(ui->lineEdit_gritSize->text().toDouble())); // [m]
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
    if (!(ui->lineEdit_diffThickness->text().isEmpty()) && ((ui->lineEdit_diffThickness->text().toDouble()) > 0)){
        return (0.001*(ui->lineEdit_diffThickness->text().toDouble())); // [m]
    } else {
        QMessageBox::warning(this, "ERROR", "Please insert a valid diffuser thickness (mm).");
        return 0;
    }
}

QString DiffuserAndObject::getGritDensity()
{
    return (ui->comboBox_gritDensity->currentText());
}
