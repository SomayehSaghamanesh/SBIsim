#include "Simulation.h"
#include "ui_Simulation.h"

#include <QMessageBox>
#include <QDebug>
#include <dlfcn.h>


Simulation::Simulation(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::Simulation)
{
    ui->setupUi(this);

    tab_source_detector = new SourceAndDetector();
    tab_diffuser_object = new DiffuserAndObject();
    tab_setup = new Setup();

    ui->SimulationTabs->addTab(tab_source_detector, QString("Source and Detector").arg(0));
    ui->SimulationTabs->addTab(tab_diffuser_object, QString("Diffuser and Object").arg(1));
    ui->SimulationTabs->addTab(tab_setup, QString("SBI Setup").arg(2));
    ui->SimulationTabs->setCurrentIndex(0);

    connect(tab_setup, &Setup::StartButtonClicked, this, &Simulation::InitiateSimulation);


}

Simulation::~Simulation()
{
    delete ui;
}


double Simulation::getFOVmin()
{

}

void Simulation::PropagateInMaterial()
{
    std::vector<std::vector<std::vector<double>>> beta, delta, t, ;

}

// sets parameters
void Simulation::setParams()
{
    // X-ray energy
    m_energyVector.clear();
    m_spectrumVector.clear();
    tab_source_detector->getXrayEnergy(m_energyVector, m_spectrumVector);

    // physics propeties of object, diffuser, and detector
    m_physList.clear();
    m_physList = mat->PhysProps(m_matList);

    // distances
    tab_setup->getDistances();
    m_SOD = tab_setup->m_SOD;
    m_SDD = tab_setup->m_SDD;
    m_SdD = tab_setup->m_SdD;


}

void Simulation::InitiateSimulation()
{
    setParams();
    // CheckUserInputs();




}
