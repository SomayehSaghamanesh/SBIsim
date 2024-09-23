#ifndef IMAGEWRITER_H
#define IMAGEWRITER_H

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include <vector>
#include <string>

template <typename T>
class ImageWriter {
public:
    using ImageType = itk::Image<T, 3>;

    ImageWriter(size_t xSize, size_t ySize, size_t zSize);
    void SetData(const std::vector<std::vector<std::vector<T>>>& data);
    void Save(const std::string& filename);

private:
    typename ImageType::Pointer m_Image;
};

#include "ImageWriter.tpp"

#endif // IMAGEWRITER_H
