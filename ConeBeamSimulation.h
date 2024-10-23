// #ifndef CONEBEAMSIMULATION_H
// #define CONEBEAMSIMULATION_H

// // #include "Compute.h"
// #include "Simulation.h"

// #include <vector>

// class ConeBeamSimulation
// {
// public:

//     ConeBeamSimulation(Simulation* sim = nullptr);

//     ~ConeBeamSimulation();

//     const double pi = M_PI;

//     void PropagateInMaterial(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const std::vector<std::vector<std::vector<int>>>& subject, const double& sliceThickness,
//                                          const double& sourceToSubjectDist, const std::vector<double>& Mag, const std::vector<double>& fresnelMag,const int& numMVSlices,
//                              const int& numVoxelSlicesInZ, const Materials::refractiveIndex& refrIndx1, const Materials::refractiveIndex& refrIndx2, const size_t& energyIndx);



//     void PropagateInFreeSpace(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const double& thickness1, const double& dist1, const double& dist2,
//                                   const double& mag2, const double& mag1, const size_t& energyIndx, const bool isLastPropagation);

//     void InitiateSimulation();

// private:

//     Simulation *simPtr;
//     // Compute *compute;
//     int m_numPixels, m_numObjVoxelsZ, m_numMVSlicesObj, m_numMVSlicesDiff, m_numDiffVoxelsZ, m_numInterp, m_numProj;
//     double m_pixelSize, m_objThickness, m_diffThickness, m_SOD, m_SdD, m_SDD;


//     void Gradient2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& Gx, std::vector<std::vector<double>>& Gy, const double spacing);

//     void Laplacian2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& L, const double& spacing);

//     std::vector<std::vector<std::vector<int>>> rotate3DObjectY(std::vector<std::vector<std::vector<int>>>& object, const int& numProjs);
//     double Interpolation2D(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y,
//                                     const std::vector<std::vector<std::vector<int>>>& grid, const int& sliceNum, const int& sliceStep, const double& n1_part,
//                                     const double& n2_part, const double& mag, const double x, const double y, double default_value);

//     void SubPixelInterpolation(const std::vector<std::vector<std::vector<int>>>& subject, std::vector<std::vector<double>>& beta,
//                                         std::vector<std::vector<double>>& delta, const Materials::refractiveIndex& n1, const Materials::refractiveIndex& n2,
//                                         const double& sourceToSubjectDist, const std::vector<double>& Mag, const int& sliceNum, const int& sliceStep, const int& m_numInterp);


// };

// #endif // CONEBEAMSIMULATION_H
