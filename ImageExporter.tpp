#include "ImageExporter.h"
#include <iostream>
#include <itkImageFileWriter.h>

template <typename T>
ImageExporter<T>::ImageExporter(size_t xSize, size_t ySize, size_t zSize) {
    typename ImageType::SizeType size;
    size[0] = xSize;
    size[1] = ySize;
    size[2] = zSize;

    typename ImageType::RegionType region;
    region.SetSize(size);

    m_Image = ImageType::New();
    m_Image->SetRegions(region);
    m_Image->Allocate();
    m_Image->FillBuffer(0);
}

template <typename T>
void ImageExporter<T>::SetData(const std::vector<std::vector<std::vector<T>>>& data) {
    if (data.size() != m_Image->GetLargestPossibleRegion().GetSize()[2] ||
        data[0].size() != m_Image->GetLargestPossibleRegion().GetSize()[1] ||
        data[0][0].size() != m_Image->GetLargestPossibleRegion().GetSize()[0]) {
        throw std::runtime_error("Data dimensions do not match image dimensions.");
    }

    typename ImageType::IndexType index;
    for (size_t z = 0; z < m_Image->GetLargestPossibleRegion().GetSize()[2]; ++z) {
        for (size_t y = 0; y < m_Image->GetLargestPossibleRegion().GetSize()[1]; ++y) {
            for (size_t x = 0; x < m_Image->GetLargestPossibleRegion().GetSize()[0]; ++x) {
                index[0] = x;
                index[1] = y;
                index[2] = z;
                m_Image->SetPixel(index, data[z][y][x]);
            }
        }
    }
}

template <typename T>
void ImageExporter<T>::Save(const std::string& filename) {
    using WriterType = itk::ImageFileWriter<ImageType>;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(m_Image);

    try {
        writer->Update();
        std::cout << "Image successfully saved as " << filename << std::endl;
    } catch (itk::ExceptionObject &error) {
        std::cerr << "Error: " << error << std::endl;
        throw;
    }
}
