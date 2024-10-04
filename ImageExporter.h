#ifndef IMAGEEXPORTER_H
#define IMAGEEXPORTER_H

#include "itkImage.h"
#include <vector>
#include <string>


template <typename T>
class ImageExporter {

public:
    using ImageType3D = itk::Image<T, 3>;
    using ImageType2D = itk::Image<T, 2>;

    ImageExporter(size_t xSize, size_t ySize, size_t zSize); // overloaded constructor
    ImageExporter(size_t xSize, size_t ySize); // overloaded constructor
    void SetData(const std::vector<std::vector<std::vector<T>>>& data); // overloaded method
    void SetData(const std::vector<std::vector<T>>& data); // overloaded method
    void Save3D(const std::string& filename);
    void Save2D(const std::string& filename);

private:
    typename ImageType3D::Pointer m_Image3D;
    typename ImageType2D::Pointer m_Image2D;

};

#include "ImageExporter.tpp"

#endif // IMAGEEXPORTER_H
