#ifndef SIMULATION_H
#define SIMULATION_H

#include "SourceAndDetector.h"
#include "DiffuserAndObject.h"
#include "Setup.h"
#include "Materials.h"

#include <QMainWindow>


namespace Ui {
class Simulation;
}

class Simulation : public QMainWindow
{
    Q_OBJECT

public:
    explicit Simulation(QWidget *parent = nullptr);
    ~Simulation();

private:
    Ui::Simulation *ui;

    SourceAndDetector *tab_source_detector;
    DiffuserAndObject *tab_diffuser_object;
    Setup *tab_setup;
    Materials *mat;

    const double hc = 1.24e-9; // [keV.m]
    const double pi = M_PI;

    std::vector<double> m_energyVector{0}, m_waveNumber{0}, m_spectrumVector{0}, m_detEnergyResponse{0};
    QVector<QString> m_matList{"", "", "", ""};
    QVector<Materials::MaterialProperties> m_physList = {{0, 0, "", ""}}; // it will have four elements at end: the order is detector, diffuser, base, and object

    double m_SOD{0}, m_SdD{0}, m_SDD{0}, m_pixelSize{0};
    int m_numProj{1}, m_numObjVoxelsZ{1}, m_numDiffVoxelsZ{1}, m_numInterp{1}, m_numMVSlicesObj{1}, m_numMVSlicesDiff{1}, m_numPixels{0};
    double m_objThickness{0}, m_diffThickness{0};
    std::vector<double> m_M_obj{1}, m_M_diff{1}, m_fM_obj{1}, m_fM_diff{1};

    // coordinates
    std::vector<std::vector<double>> m_X, m_Y, m_rsqr;

    // Images
    std::vector<std::vector<double>> m_I_bg, m_I_fg, m_phi_bg, m_phi_fg;

    // tempoarary loop matrices
    std::vector<std::vector<double>> beta; // imaginary part of the complex refractive index
    std::vector<std::vector<double>> delta; // real part of the complex refractive index
    std::vector<std::vector<double>> Gx_I; // intensity gradient along x : here vertical side
    std::vector<std::vector<double>> Gy_I; // intensity gradient along y : here horizontal side
    std::vector<std::vector<double>> Gx_phi; // phase gradient along x : here vertical side
    std::vector<std::vector<double>> Gy_phi; // phase gradient along y : here horizontal side
    std::vector<std::vector<double>> L_phi; // 2D Laplacian of phi

    void getDynamicMagFactors();
    void setParams();

    std::vector<std::vector<std::vector<int>>> rotate3DObjectY(std::vector<std::vector<std::vector<int>>>& object, const int& numProjs);

    void Gradient2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& Gx, std::vector<std::vector<double>>& Gy, const double spacing);
    void Laplacian2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& L, const double& spacing);

    double Interpolation2D(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y,
                           const std::vector<std::vector<std::vector<int>>>& grid, const int& sliceNum, const int& sliceStep, const double& n1_part,
                           const double& n2_part, const double& mag, const double x, const double y, double default_value = 0);

    void SubPixelInterpolation(const std::vector<std::vector<std::vector<int>>>& subject, std::vector<std::vector<double>>& beta,
                               std::vector<std::vector<double>>& delta, const Materials::refractiveIndex& n1, const Materials::refractiveIndex& n2,
                               const double& sourceToSubjectDist, const std::vector<double>& Mag, const int& sliceNum, const int& sliceStep, const int& m_numInterp);

    void PropagateInMaterial(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const std::vector<std::vector<std::vector<int>>>& subject,
                             const double& thickness, const double& sourceToSubjectDist, const std::vector<double>& Mag,
                             const std::vector<double>& fresnelMag,const int& numMVSlices, const int& numVoxelSlicesInZ,
                             const Materials::refractiveIndex& refrIndx1, const Materials::refractiveIndex& refrIndx2, const size_t& energyIndx);

    void PropagateInFreeSpace(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const double& thickness1,
                              const double& dist1, const double& dist2, const double& Mag1,
                              const double& Mag2, const size_t& energyIndx, const bool isLastPropagation);
    void ConeBeamSimulation();
    void ParBeamSimulation();
    void InitiateSimulation();

};

#endif // SIMULATION_H
