#ifndef DIFFUSER_H
#define DIFFUSER_H

#include "DiffuserAndObject.h"
// #include "SourceAndDetector.h"
#include "Setup.h"
#include <random>

class Diffuser : public DiffuserAndObject
{

public:

    Diffuser();
    ~Diffuser();

    void CreateDiffuser(std::vector<std::vector<std::vector<double>>>& diffuser, double refrIndx);

private:

    std::array<int, 3> m_diffuserSize;
    int m_numGrit;
    void getRandomValue(std::vector<int>& myVec, int lower_indx, int upper_indx);
    void DistributeGrits(std::vector<int>& cx, std::vector<int>& cy, std::vector<int>& cz, std::vector<int>& r);


};

#endif // DIFFUSER_H
