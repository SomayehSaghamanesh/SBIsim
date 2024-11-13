#ifndef OBJECT_H
#define OBJECT_H

#include "DiffuserAndObject.h"
#include <cmath>

class Object
{
public:
    Object(int numObjVoxelsZ, int numPixels, double pixelObj, DiffuserAndObject* diffuserAndObject);
    ~Object();

    void CreateObject(std::vector<std::vector<std::vector<int>>>& object);
    void CreateObject(std::vector<std::vector<std::vector<float>>>& object); // overloaded
    std::vector<std::vector<std::vector<int>>> rotate3DObjectZ(std::vector<std::vector<std::vector<int>>>& object);

private:

    int m_numObjVoxelsZ, m_numPixels;
    double m_pixelObj;
    DiffuserAndObject* m_diffuserAndObject;

    void CreateSphere(std::vector<std::vector<std::vector<int>>>& object);
    void CreateCylinder(std::vector<std::vector<std::vector<int>>>& object);

};

#endif // OBJECT_H
