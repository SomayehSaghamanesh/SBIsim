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

    std::vector<double> m_energyVector{0};
    std::vector<double> m_spectrumVector{0};
    QVector<QString> m_matList{"", "", "", ""};
    QVector<Materials::MaterialProperties> m_physList = {{0, 0, "", ""}}; // it will have four elements at end: the order is detector, diffuser, base, and object

    double m_SOD{0}, m_SdD{0}, m_SDD{0}, m_pixelSize{0};
    int m_numProj{1}, m_numObjVoxelsZ{1}, m_numDiffVoxelsZ{1}, m_numInterp{1}, m_numMVSlicesObj{1}, m_numMVSlicesDiff{1}, m_numPixels{0}, m_objectType{0};
    double m_M_obj{1}, m_M_diff{1}, m_objThickness{0}, m_diffThickness{0}, m_pixelObj{0}, m_pixelDiff{0};
    std::vector<double> m_dM_obj{1}, m_dM_diff{1}, m_fM_obj{1}, m_fM_diff{1};

    // coordinates
    std::vector<std::vector<float>> m_X, m_Y, m_rsqr;

    // Images
    std::vector<std::vector<float>> m_I, m_phi;

    void getDynamicMagFactors();
    void setParams();
    double Interpolation2D(const std::vector<std::vector<float>>& X, const std::vector<std::vector<float>>& Y,
                           const std::vector<std::vector<std::vector<int>>>& grid, const int& sliceNum, const double& n1_part,
                           const double& n2_part, const double x, const double y, double default_value = 0);

    void SubPixelInterpolation(const std::vector<std::vector<std::vector<int>>>& subject, std::vector<std::vector<double>>& beta,
                               std::vector<std::vector<double>>& delta, const Materials::refractiveIndex& n1, const Materials::refractiveIndex& n2,
                               const double& thickness, const double& distToDet, const double& globalMag, const int& sliceNum, const int& m_numInterp);

    void PropagateInMaterial(const std::vector<std::vector<std::vector<int>>>& subject, const double& thickness, const double& distToDet, const double& globalMag, const std::vector<double>& dynamicMag,
                             const int& numMVSlices, const int& numVoxelSlicesInZ, const Materials::refractiveIndex& refrIndx1, const Materials::refractiveIndex& refrIndx2);
    void PropagateInFreeSpace();
    void InitiateSimulation();

};

#endif // SIMULATION_H
