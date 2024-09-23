#include "ImageWriter.h"
#include "itkImageFileWriter.h"
#include "itkExceptionObject.h"
#include <iostream>

// Constructor: Initializes image with given dimensions
ImageWriter::ImageWriter(size_t xSize, size_t ySize, size_t zSize)
    : xSize(xSize), ySize(ySize), zSize(zSize) {
    InitializeImage();
}

// Private method to initialize the image with the specified size
void ImageWriter::InitializeImage() {
    image = ImageType::New();
    ImageType::SizeType size;
    size[0] = xSize;
    size[1] = ySize;
    size[2] = zSize;

    ImageType::RegionType region;
    region.SetSize(size);

    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0); // Optional: initialize all pixels to zero
}

// Set image data from a provided 1D vector
void ImageWriter::SetImageData(const std::vector<float>& data) {
    if (data.size() != xSize * ySize * zSize) {
        std::cerr << "Error: Data size does not match image dimensions!" << std::endl;
        return;
    }

    ImageType::IndexType index;
    size_t counter = 0;
    for (size_t z = 0; z < zSize; ++z) {
        for (size_t y = 0; y < ySize; ++y) {
            for (size_t x = 0; x < xSize; ++x) {
                index[0] = x;
                index[1] = y;
                index[2] = z;
                image->SetPixel(index, data[counter++]);
            }
        }
    }
}

// Save the image to the specified filename; format inferred from file extension
bool ImageWriter::SaveImage(const std::string& filename) {
    using WriterType = itk::ImageFileWriter<ImageType>;
    WriterType::Pointer writer = WriterType::New();

    writer->SetFileName(filename); // Set the output filename
    writer->SetInput(image);       // Connect the image to the writer

    try {
        writer->Update(); // Write the image to the file
        std::cout << "Image successfully saved as " << filename << std::endl;
        return true;
    } catch (itk::ExceptionObject& error) {
        std::cerr << "Error saving image: " << error << std::endl;
        return false;
    }
}

