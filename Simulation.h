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

    // bool m_radioFlag{false}, m_tomoFlag{false}, m_parBeamFlag{false}, m_coneBeamFlag{false}, m_monochrom{false}, m_polychromFlag{false};

    const double hc = 1.24e-9; // [keV.m]
    const double pi = M_PI;

    std::vector<double> m_energyVector{0};
    std::vector<double> m_spectrumVector{0};
    QVector<QString> m_matList{"", "", "", ""};
    // QVector<QString> m_matList{tab_source_detector->getDetectorMaterial(), tab_diffuser_object->getGritMaterial(),
                               // tab_diffuser_object->getBaseMaterial(), tab_diffuser_object->getObjectMaterial()};
    QVector<Materials::MaterialProperties> m_physList = {{0, 0, "", ""}}; // it will have four elements at end: the order is detector, diffuser, base, and object
    int m_numProj{1};

    double m_SOD{0}, m_SdD{0}, m_SDD{0};
    double m_pixelSize{0};
    int m_numObjVoxelsZ{1}, m_numDiffVoxelsZ{1}, m_numInterp{1}, m_numMVSlicesObj{1}, m_numMVSlicesDiff{1}, m_numPixels{0};

    double m_M_obj{1}, m_M_diff{1}, m_objThickness{0}, m_diffThickness{0}, m_pixelObj{0}, m_pixelDiff{0};
    std::vector<double> m_dM_obj{1}, m_dM_diff{1}, m_fM_obj{1}, m_fM_diff{1};

    void getDynamicMagFactors();
    void setParams();
    void InitiateSimulation();
    void PropagateInMaterial();
    void PropagateInFreeSpace();

    double getFOVmin();

};

#endif // SIMULATION_H
