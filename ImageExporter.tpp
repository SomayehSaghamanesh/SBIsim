#include "ImageExporter.h"
#include <iostream>
#include <itkImageFileWriter.h>

template <typename T>
ImageExporter<T>::ImageExporter(size_t xSize, size_t ySize, size_t zSize) {

    typename ImageType3D::SizeType size;
    size[0] = xSize;
    size[1] = ySize;
    size[2] = zSize;

    typename ImageType3D::RegionType region;
    region.SetSize(size);

    m_Image3D = ImageType3D::New();
    m_Image3D->SetRegions(region);
    m_Image3D->Allocate();
    m_Image3D->FillBuffer(0);
}

template <typename T>
ImageExporter<T>::ImageExporter(size_t xSize, size_t ySize) {

    typename ImageType2D::SizeType size;
    size[0] = xSize;
    size[1] = ySize;

    typename ImageType2D::RegionType region;
    region.SetSize(size);

    m_Image2D = ImageType2D::New();
    m_Image2D->SetRegions(region);
    m_Image2D->Allocate();
    m_Image2D->FillBuffer(0);
}

template <typename T>
void ImageExporter<T>::SetData(const std::vector<std::vector<std::vector<T>>>& data) {
    if (data.size() != m_Image3D->GetLargestPossibleRegion().GetSize()[2] ||
        data[0].size() != m_Image3D->GetLargestPossibleRegion().GetSize()[1] ||
        data[0][0].size() != m_Image3D->GetLargestPossibleRegion().GetSize()[0]) {
        throw std::runtime_error("Data dimensions do not match image dimensions.");
    }

    typename ImageType3D::IndexType index;
    for (size_t z = 0; z < m_Image3D->GetLargestPossibleRegion().GetSize()[2]; ++z) {
        for (size_t y = 0; y < m_Image3D->GetLargestPossibleRegion().GetSize()[1]; ++y) {
            for (size_t x = 0; x < m_Image3D->GetLargestPossibleRegion().GetSize()[0]; ++x) {
                index[0] = x;
                index[1] = y;
                index[2] = z;
                m_Image3D->SetPixel(index, data[z][y][x]);
            }
        }
    }
}

template <typename T>
void ImageExporter<T>::SetData(const std::vector<std::vector<T>>& data) {
    if (data.size() != m_Image2D->GetLargestPossibleRegion().GetSize()[1] ||
        data[0].size() != m_Image2D->GetLargestPossibleRegion().GetSize()[0]) {
        throw std::runtime_error("Data dimensions do not match image dimensions.");
    }

    typename ImageType2D::IndexType index;
    for (size_t y = 0; y < m_Image2D->GetLargestPossibleRegion().GetSize()[1]; ++y) {
        for (size_t x = 0; x < m_Image2D->GetLargestPossibleRegion().GetSize()[0]; ++x) {
            index[0] = x;
            index[1] = y;
            m_Image2D->SetPixel(index, data[y][x]);
        }
    }
}

template <typename T>
void ImageExporter<T>::Save3D(const std::string& filename) {
    using WriterType = itk::ImageFileWriter<ImageType3D>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(m_Image3D);

    try {
        writer->Update();
        std::cout << "Image successfully saved as " << filename << std::endl;
    } catch (itk::ExceptionObject &error) {
        std::cerr << "Error: " << error << std::endl;
        throw;
    }
}

template <typename T>
void ImageExporter<T>::Save2D(const std::string& filename) {
    using WriterType = itk::ImageFileWriter<ImageType2D>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(m_Image2D);

    try {
        writer->Update();
        std::cout << "Image successfully saved as " << filename << std::endl;
    } catch (itk::ExceptionObject &error) {
        std::cerr << "Error: " << error << std::endl;
        throw;
    }
}
