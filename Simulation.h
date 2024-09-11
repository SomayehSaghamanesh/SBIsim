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

    Setup *tab_setup;
    DiffuserAndObject *tab_diffuser_object;
    SourceAndDetector *tab_source_detector;
    Materials *mat;

    bool m_radioFlag{false}, m_tomoFlag{false}, m_parBeamFlag{false}, m_coneBeamFlag{false}, m_monochrom{false}, m_polychromFlag{false};

    const double hc = 1.24e-9; // [keV.m]
    const double pi = M_PI;

    std::vector<double> m_energyVector{0};
    std::vector<double> m_spectrumVector{0};
    QString m_objectMat{""}, m_gritMat{""}, m_baseMat{""}, m_detectorMat{""};
    QVector<QString> m_matList{m_objectMat, m_gritMat, m_baseMat, m_detectorMat};
    QVector<Materials::MaterialProperties> m_physList = {{0, 0, "", ""}};


    double m_SOD{0}, m_SdD{0}, m_SDD{0};

    void setParams();
    void InitiateSimulation();
    void PropagateInMaterial();
    void PropagateInFreeSpace();

    double getFOVmin();

};

#endif // SIMULATION_H
