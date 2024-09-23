#include "Simulation.h"
#include "ui_Simulation.h"
#include "Object.h"
#include "Diffuser.h"

#include <QMessageBox>
#include <QDebug>
#include <dlfcn.h>
#include <memory>
#include <fstream>


Simulation::Simulation(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::Simulation)
{
    ui->setupUi(this);

    tab_source_detector = new SourceAndDetector();
    tab_diffuser_object = new DiffuserAndObject();
    tab_setup = new Setup();
    mat = new Materials();

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

// double Simulation::getFOVmin()
// {

// }

void Simulation::getDynamicMagFactors()
{
    m_fM_obj.clear();
    m_fM_diff.clear();

    // object and diff magnification factors at each slice inside obj/diff
    m_dM_obj.clear();
    m_dM_diff.clear();

    for (int i = 0 ; i <m_numObjVoxelsZ ; i++)
    {
        m_dM_obj.push_back(m_SDD/(m_SOD - m_objThickness/2 + (i - 1)*m_pixelSize/m_M_obj)); // dynamic magnification factors of slices
        m_fM_obj.push_back((m_SOD - m_objThickness/2 + (i - 1)*m_pixelSize/m_M_obj)/(m_SOD - m_objThickness/2));; // Magnification after applying Fresnel scaling theory
    }

    for (int i = 0 ; i <m_numDiffVoxelsZ ; i++)
    {
        m_dM_diff.push_back(m_SDD/(m_SdD - m_diffThickness/2 + (i - 1)*m_pixelSize/m_M_diff)); // dynamic magnification factors of slices
        m_fM_diff.push_back((m_SdD - m_diffThickness/2 + (i - 1)*m_pixelSize/m_M_diff)/(m_SdD - m_diffThickness/2)); // Magnification after applying Fresnel scaling theory
    }
}

void Simulation::PropagateInMaterial()
{
    std::vector<std::vector<std::vector<double>>> beta, delta, t;
}

// sets parameters
void Simulation::setParams()
{
    // tomo
    m_numProj = tab_setup->getNumProjections();

    // X-ray energy
    m_energyVector.clear();
    m_spectrumVector.clear();
    tab_source_detector->getXrayEnergy(m_energyVector, m_spectrumVector);

    // detector
    tab_source_detector->getPixelSizeMM();
    m_pixelSize = tab_source_detector->m_pixelSizeMM;
    m_numPixels = tab_source_detector->getNumPixels();
    //det thickness?

    // physics propeties of object, diffuser, and detector
    m_matList.clear();
    m_matList = {tab_source_detector->getDetectorMaterial(), tab_diffuser_object->getGritMaterial(),
                tab_diffuser_object->getBaseMaterial(), tab_diffuser_object->getObjectMaterial()};
    m_physList.clear();
    m_physList = mat->PhysProps(m_matList);

    // distances
    tab_setup->getDistances();
    m_SOD = tab_setup->m_SOD;
    m_SdD = tab_setup->m_SdD;
    m_SDD = tab_setup->m_SDD;

    // global magnification factors
    m_M_obj = m_SDD/m_SOD;
    m_M_diff = m_SDD/m_SdD;
    m_pixelObj = m_pixelSize/m_M_obj;
    m_pixelDiff = m_pixelSize/m_M_diff;

    m_objThickness = tab_diffuser_object->getObjectThickness();
    m_diffThickness = tab_diffuser_object->getDiffuserThickness();

    m_numMVSlicesObj = tab_diffuser_object->getNumMVSlices(m_objThickness, m_numObjVoxelsZ, m_pixelObj);
    m_numMVSlicesDiff = tab_diffuser_object->getNumMVSlices(m_diffThickness, m_numDiffVoxelsZ, m_pixelDiff);
    m_numInterp = tab_diffuser_object->getNumInterpolations();

    // mag factors based on the fresnel scaling theory
    getDynamicMagFactors();

    // for (uint i = 0 ; i < m_energyVector.size() ; i++)
    // {
    //     qDebug() << "E : " << m_energyVector[i] << "\t" << m_spectrumVector[i];
    // }

    qDebug() << "pixel size : " << m_pixelSize;
    qDebug() << "SOD : " << m_SOD;
    qDebug() << "SdD : " << m_SdD;
    qDebug() << "SDD : " << m_SDD;

    for (uint i = 0 ; i < m_physList.size() ; i++)
    {
        qDebug() << "physlist : d = " << m_physList[i].density << ", Z/A = " << m_physList[i].Z_A << ", name = " << m_physList[i].name;
    }

    qDebug() << "numPixels : " << m_numPixels;
    qDebug() << "obj-thickness : " << m_objThickness;
    qDebug() << "numVoxelsInZObj : " << m_numObjVoxelsZ;
    qDebug() << "numMVSlices_obj : " << m_numMVSlicesObj;
    qDebug() << "numMVSlices_diff : " << m_numMVSlicesDiff;
    qDebug() << "numInterp : " << m_numInterp;
    qDebug() << "numProj : " << m_numProj;
    qDebug() << "M_obj : " << m_M_obj;
    qDebug() << "M_diff : " << m_M_diff;
    qDebug() << "MM_obj : " << m_dM_obj;
    qDebug() << "MM_diff : " << m_dM_diff;
    qDebug() << "m_pixelObj : " << m_pixelObj;
    qDebug() << "m_pixelDiff : " << m_pixelDiff;
}

void Simulation::InitiateSimulation()
{
    setParams();
    // CheckUserInputs();

    // Create object
    std::vector<std::vector<std::vector<Materials::refractiveIndex>>> obj(m_numPixels, std::vector<std::vector<Materials::refractiveIndex>>(m_numPixels, std::vector<Materials::refractiveIndex>(m_numObjVoxelsZ, Materials::refractiveIndex(0.0, 0.0))));
    std::vector<Materials::refractiveIndex> n_obj = mat->RefractiveIndex(m_energyVector, m_physList[3].formula, m_physList[3].density);
    std::unique_ptr<Object> object = std::make_unique<Object>(m_numObjVoxelsZ, m_numPixels, m_pixelObj, tab_diffuser_object);

    // Create diffuser
    std::vector<std::vector<std::vector<Materials::refractiveIndex>>> diff(m_numPixels, std::vector<std::vector<Materials::refractiveIndex>>(m_numPixels, std::vector<Materials::refractiveIndex>(m_numDiffVoxelsZ, Materials::refractiveIndex(0.0, 0.0))));
    std::vector<Materials::refractiveIndex> n_diff = mat->RefractiveIndex(m_energyVector, m_physList[1].formula, m_physList[1].density);
    std::vector<Materials::refractiveIndex> n_base = mat->RefractiveIndex(m_energyVector, m_physList[2].formula, m_physList[2].density);
    std::unique_ptr<Diffuser> diffuser = std::make_unique<Diffuser>(m_numPixels, m_numDiffVoxelsZ, m_pixelDiff, tab_diffuser_object);

    // loop over projections
    for (int proj = 0 ; proj < m_numProj ; proj++)
    {
        // loop on all X-ray energies
        for (size_t i = 0 ; i < m_energyVector.size() ; i++)
        {
            if (tab_diffuser_object->m_sphere){
                object->CreateSphere(obj, n_obj[i]);

            } else if (tab_diffuser_object->m_cylinder){
                object->CreateCylinder(obj, n_obj[i]);

            } else if (tab_diffuser_object->m_virtualObject){
                //m_vObject
            }

            diffuser->CreateDiffuser(diff, n_diff[i], n_base[i]);

            // std::ofstream file_obj("object.txt");
            // std::ofstream file_diff("diffuser.txt");

            // if ((file_obj.is_open()) && (file_diff.is_open())){
            //     for (int p = 0 ; p < m_numPixels ; p++)
            //     {
            //         for (int q = 0 ; q < m_numPixels ; q++)
            //         {
            //             for (int l = 0 ; l < m_numObjVoxelsZ; l++)
            //             {
            //                 file_obj << obj[p][q][l].imaginaryPart << " ";
            //             }
            //             file_obj << "\n";
            //             for (int l = 0 ; l < m_numDiffVoxelsZ; l++)
            //             {
            //                 file_diff << diff[p][q][l].imaginaryPart << " ";
            //             }
            //             file_diff << "\n";
            //         }
            //     }
            // } else {
            //     qDebug() << "file was not open";
            //     return;
            // }

            // file_obj.close();
            // file_diff.close();

            // qDebug() << "n_obj-im : " << n_obj[i].imaginaryPart << "; n_obj-re : " << n_obj[i].realPart;
            // qDebug() << "n_diff-im : " << n_diff[i].imaginaryPart << "; n_diff-re : " << n_diff[i].realPart;

            // must reset obj/diff matrices to 0 to recreate them
            for (int m = 0 ; m < m_numPixels ; m++)
            {
                for (int n = 0 ; n < m_numPixels ; n++)
                {
                    obj[m][n].assign(m_numObjVoxelsZ, Materials::refractiveIndex(0.0, 0.0));
                    diff[m][n].assign(m_numDiffVoxelsZ, Materials::refractiveIndex(0.0, 0.0));
                }
            }
        }
    }
}
