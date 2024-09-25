#ifndef DIFFUSER_H
#define DIFFUSER_H

#include "DiffuserAndObject.h"

class Diffuser
{

public:

    Diffuser(int numPixels, int numDiffVoxelsZ, double pixelDiff, DiffuserAndObject *diffuserAndObject);
    ~Diffuser();

    void CreateDiffuser(std::vector<std::vector<std::vector<int>>>& diffuser);

private:

    int m_numPixels, m_numDiffVoxelsZ;
    double m_pixelDiff;
    DiffuserAndObject* m_diffuserAndObject;
    std::array<int, 3> m_diffuserSize;
    int m_numGrit;

    void getRandomValue(std::vector<int>& myVec, int lower_indx, int upper_indx);
    void DistributeGrits(std::vector<int>& cx, std::vector<int>& cy, std::vector<int>& cz, std::vector<int>& r);

};

#endif // DIFFUSER_H
