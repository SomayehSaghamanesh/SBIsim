#ifndef IMAGEEXPORTER_H
#define IMAGEEXPORTER_H

#include "itkImage.h"
#include <vector>
#include <string>


template <typename T>
class ImageExporter {

public:
    using ImageType = itk::Image<T, 3>;

    ImageExporter(size_t xSize, size_t ySize, size_t zSize);
    void SetData(const std::vector<std::vector<std::vector<T>>>& data);
    void Save(const std::string& filename);

private:
    typename ImageType::Pointer m_Image;

};

// #include "ImageExporter.tpp"


#endif // IMAGEEXPORTER_H
