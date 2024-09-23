#ifndef OBJECT_H
#define OBJECT_H

#include "DiffuserAndObject.h"
#include "Materials.h"

class Object
{
public:
    Object(int numObjVoxelsZ, int numPixels, double pixelObj, DiffuserAndObject* diffuserAndObject);
    ~Object();

    void CreateSphere(std::vector<std::vector<std::vector<Materials::refractiveIndex>>>& object, Materials::refractiveIndex& n_obj);
    void CreateCylinder(std::vector<std::vector<std::vector<Materials::refractiveIndex>>>& object, Materials::refractiveIndex& n_obj);

private:

    int m_numObjVoxelsZ, m_numPixels;
    double m_pixelObj;
    DiffuserAndObject* m_diffuserAndObject;


};

#endif // OBJECT_H
