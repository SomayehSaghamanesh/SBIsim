#ifndef OBJECT_H
#define OBJECT_H

#include "DiffuserAndObject.h"
#include "SourceAndDetector.h"
#include "Setup.h"

class Object : public DiffuserAndObject, public SourceAndDetector
{
public:
    Object();
    ~Object();

    void CreateSphere(std::vector<std::vector<std::vector<double>>>& object, double refrIndx);
    void CreateCylinder(std::vector<std::vector<std::vector<double>>>& object, double refrIndx);

private:


};

#endif // OBJECT_H
