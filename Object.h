#ifndef OBJECT_H
#define OBJECT_H

#include "DiffuserAndObject.h"

class Object
{
public:
    Object(int numObjVoxelsZ, int numPixels, double pixelObj, DiffuserAndObject* diffuserAndObject, int objectType);
    ~Object();

    void CreateObject(std::vector<std::vector<std::vector<int>>>& object);

private:

    int m_numObjVoxelsZ, m_numPixels;
    double m_pixelObj;
    DiffuserAndObject* m_diffuserAndObject;
    int m_objectType;

    void CreateSphere(std::vector<std::vector<std::vector<int>>>& object);
    void CreateCylinder(std::vector<std::vector<std::vector<int>>>& object);
    void CreateVirtualObject();


};

#endif // OBJECT_H
