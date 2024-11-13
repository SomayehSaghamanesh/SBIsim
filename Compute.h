#ifndef COMPUTE_H
#define COMPUTE_H

#include "Materials.h"

#include <vector>
#include <cmath>


class Compute
{
public:
    Compute();
    ~Compute();

    const double pi = M_PI;

    void Gradient2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& Gx, std::vector<std::vector<double>>& Gy, const double spacing);


    void Laplacian2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& L, const double& spacing);


    std::vector<std::vector<std::vector<int>>> rotate3DObjectY(std::vector<std::vector<std::vector<int>>>& object, const int& numProjs);


    double Interpolation2D(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y,
                                    const std::vector<std::vector<std::vector<int>>>& grid, const int& sliceNum, const int& sliceStep, const double& n1_part,
                                    const double& n2_part, const double& mag, const double x, const double y, double default_value);


    void SubPixelInterpolation(const std::vector<std::vector<std::vector<int>>>& subject, std::vector<std::vector<double>>& beta,
                                        std::vector<std::vector<double>>& delta, const Materials::refractiveIndex& n1, const Materials::refractiveIndex& n2,
                               const double& sourceToSubjectDist, const std::vector<double>& Mag, const int& sliceNum, const int& sliceStep, const int& m_numInterp);

};

#endif // COMPUTE_H
