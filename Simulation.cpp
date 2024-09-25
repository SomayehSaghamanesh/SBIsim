#include "Simulation.h"
#include "ui_Simulation.h"
#include "Object.h"
#include "Diffuser.h"

#include <QMessageBox>
#include <QDebug>
#include <dlfcn.h>
#include <memory>
// #include <iostream>


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

void Simulation::setParams()
{
    // tomo
    m_numProj = tab_setup->getNumProjections();

    // X-ray energy
    m_energyVector.clear();
    m_spectrumVector.clear();
    tab_source_detector->getXrayEnergy(m_energyVector, m_spectrumVector);

    // detector
    m_pixelSize = tab_source_detector->getPixelSizeMM();
    m_numPixels = tab_source_detector->getNumPixels();
    tab_source_detector->Meshgrid(m_X, m_Y, m_numPixels, m_pixelSize);
    tab_source_detector->DetectorCoordinates(m_X, m_Y, m_rsqr, m_numPixels);
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

    // object type
    if (tab_diffuser_object->m_sphere){
        m_objectType = 1;

    } else if (tab_diffuser_object->m_cylinder){
        m_objectType = 2;

    } else if (tab_diffuser_object->m_virtualObject){
        m_objectType = 3;
    }

    // below part only for debugging
    // for (uint i = 0 ; i < m_energyVector.size() ; i++)
    // {
    //     qDebug() << "E : " << m_energyVector[i] << "\t" << m_spectrumVector[i];
    // }

    // qDebug() << "pixel size : " << m_pixelSize;
    // qDebug() << "SOD : " << m_SOD;
    // qDebug() << "SdD : " << m_SdD;
    // qDebug() << "SDD : " << m_SDD;

    // for (uint i = 0 ; i < m_physList.size() ; i++)
    // {
    //     qDebug() << "physlist : d = " << m_physList[i].density << ", Z/A = " << m_physList[i].Z_A << ", name = " << m_physList[i].name;
    // }

    // qDebug() << "numPixels : " << m_numPixels;
    // qDebug() << "obj-thickness : " << m_objThickness;
    // qDebug() << "numVoxelsInZObj : " << m_numObjVoxelsZ;
    // qDebug() << "numMVSlices_obj : " << m_numMVSlicesObj;
    // qDebug() << "numMVSlices_diff : " << m_numMVSlicesDiff;
    // qDebug() << "numInterp : " << m_numInterp;
    // qDebug() << "numProj : " << m_numProj;
    // qDebug() << "M_obj : " << m_M_obj;
    // qDebug() << "M_diff : " << m_M_diff;
    // qDebug() << "MM_obj : " << m_dM_obj;
    // qDebug() << "MM_diff : " << m_dM_diff;
    // qDebug() << "m_pixelObj : " << m_pixelObj;
    // qDebug() << "m_pixelDiff : " << m_pixelDiff;
}

double Simulation::Interpolation2D(const std::vector<std::vector<float>>& X, const std::vector<std::vector<float>>& Y,
                                   const std::vector<std::vector<std::vector<int>>>& grid, const int& sliceNum, const double& n1_part,
                                   const double& n2_part, const double x, const double y, double default_value)
{
    int rows = grid.size();
    int cols = grid[0].size();

    // Find the four surrounding points in X and Y
    int x_low = 0, x_high = 0, y_low = 0, y_high = 0;

    // Find the correct x indices (x_low and x_high)
    for (int i = 0; i < cols - 1; ++i) {
        if (X[i][0] <= x && X[i + 1][0] >= x) {
            x_low = i;
            x_high = i + 1;
            break;
        }
    }
    // Find the correct y indices (y_low and y_high)
    for (int i = 0; i < rows - 1; ++i) {
        if (Y[0][i] <= y && Y[0][i + 1] >= y) {
            y_low = i;
            y_high = i + 1;
            break;
        }
    }
    // Check bounds, if out of range, return default value
    if (x_low < 0 || x_high >= rows || y_low < 0 || y_high >= cols || x_low == x_high || y_low == y_high) {
        return default_value;
    }

    // Bilinear interpolation
    double x_ratio = (x - X[x_low][0]) / (X[x_high][0] - X[x_low][0]);
    double y_ratio = (y - Y[0][y_low]) / (Y[0][y_high] - Y[0][y_low]);

    double top = ( (grid[x_low][y_low][sliceNum] == 1) ? ((n1_part) * grid[x_low][y_low][sliceNum]) : ((n2_part/2) * grid[x_low][y_low][sliceNum]) ) * (1 - y_ratio)
                 + ( (grid[x_low][y_high][sliceNum] == 1) ? ((n1_part) * grid[x_low][y_high][sliceNum]) : ((n2_part/2) * grid[x_low][y_high][sliceNum]) ) * y_ratio;

    double bottom = ( (grid[x_high][y_low][sliceNum] == 1) ? ((n1_part) * grid[x_high][y_low][sliceNum]) : ((n2_part/2) * grid[x_high][y_low][sliceNum]) ) * (1 - y_ratio)
                    + ( (grid[x_high][y_high][sliceNum] == 1) ? ((n1_part) * grid[x_high][y_high][sliceNum]) : ((n2_part/2) * grid[x_high][y_high][sliceNum]) ) * y_ratio;

        return (1 - x_ratio) * top + x_ratio * bottom;

}

void Simulation::SubPixelInterpolation(const std::vector<std::vector<std::vector<int>>>& subject, std::vector<std::vector<double>>& beta,
                                       std::vector<std::vector<double>>& delta, const Materials::refractiveIndex& n1, const Materials::refractiveIndex& n2,
                                       const double& thickness, const double& distToDet, const double& globalMag, const int& sliceNum, const int& m_numInterp)
{
    double M_interp = 0;
    // loop over different interpolation steps and accumulating them
    for (int interp = 0 ; interp < m_numInterp ; interp++)
    {
        M_interp = m_SDD/(distToDet - thickness/2 + (sliceNum-1)*m_pixelSize/globalMag + interp*m_pixelSize/globalMag/m_numInterp);

        for (int i = 0 ; i < m_numPixels ; i++){
            for (int j = 0 ; j < m_numPixels ; j++){
                beta[i][j]+=Interpolation2D(m_X, m_Y, subject, n1.imaginaryPart, n2.imaginaryPart, sliceNum, m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
                delta[i][j]+=Interpolation2D(m_X, m_Y, subject, n1.realPart, n2.realPart, sliceNum, m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
            }
        }
    }
    // averaging over interpolations
    for (int i = 0 ; i < m_numPixels ; i++){
        for (int j = 0 ; j < m_numPixels ; j++){
            beta[i][j] = beta[i][j] / m_numInterp;
            delta[i][j]= delta[i][j] / m_numInterp;
        }
    }
}

void Simulation::PropagateInFreeSpace()
{

}

void Simulation::PropagateInMaterial(const std::vector<std::vector<std::vector<int>>>& subject, const double& thickness, const double& distToDet, const double& globalMag, const std::vector<double>& dynamicMag,
                                     const int& numMVSlices, const int& numVoxelSlicesInZ, const Materials::refractiveIndex& refrIndx1, const Materials::refractiveIndex& refrIndx2)
{
    // cone-beam slice thickness at plane
    std::vector<std::vector<float>> coneBeamSliceThickness(m_numPixels, std::vector<float>(m_numPixels, 0.000000));
    std::vector<std::vector<double>> beta(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
    std::vector<std::vector<double>> delta(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));

    // create magnified thickness mesh and interpolate refrative index in sub-pixel resolution
    for (int s = 1 ; s < (numVoxelSlicesInZ - 1) ; s+=numVoxelSlicesInZ/numMVSlices)
    {
        for (int i = 0 ; i < m_numPixels ; i++)
        {
            for (int j = 0 ; j < m_numPixels ; j++)
            {
                coneBeamSliceThickness[i][j] =
                    std::sqrt( m_rsqr[i][j]/(std::pow(dynamicMag[s], 2)) + std::pow((distToDet - thickness/2 + s*m_pixelSize/globalMag), 2) )
                        - std::sqrt( m_rsqr[i][j]/(std::pow(dynamicMag[s-1], 2)) + std::pow((distToDet - thickness/2 + (s-1)*m_pixelSize/globalMag), 2) );

                SubPixelInterpolation(subject, beta, delta, refrIndx1, refrIndx2, thickness, distToDet, globalMag, s, m_numInterp);
            }
        }
        beta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
        delta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
    }
}

void Simulation::InitiateSimulation()
{
    setParams();
    // CheckUserInputs();

    // Create object
    std::vector<std::vector<std::vector<int>>> obj(m_numPixels, std::vector<std::vector<int>>(m_numPixels, std::vector<int>(m_numObjVoxelsZ, 0)));
    std::vector<Materials::refractiveIndex> n_obj = mat->RefractiveIndex(m_energyVector, m_physList[3].formula, m_physList[3].density);
    std::unique_ptr<Object> object = std::make_unique<Object>(m_numObjVoxelsZ, m_numPixels, m_pixelObj, tab_diffuser_object, m_objectType);

    // Create diffuser
    std::vector<std::vector<std::vector<int>>> diff(m_numPixels, std::vector<std::vector<int>>(m_numPixels, std::vector<int>(m_numDiffVoxelsZ, 0)));
    std::vector<Materials::refractiveIndex> n_diff = mat->RefractiveIndex(m_energyVector, m_physList[1].formula, m_physList[1].density);
    std::vector<Materials::refractiveIndex> n_base = mat->RefractiveIndex(m_energyVector, m_physList[2].formula, m_physList[2].density);
    std::unique_ptr<Diffuser> diffuser = std::make_unique<Diffuser>(m_numPixels, m_numDiffVoxelsZ, m_pixelDiff, tab_diffuser_object);

    // set intensity distribution's amplitude and phase
    std::vector<std::vector<float>> I(m_numPixels, std::vector<float>(m_numPixels, 1.0));
    std::vector<std::vector<float>> phi(m_numPixels, std::vector<float>(m_numPixels, 0.0));

    // loop over projections
    for (int proj = 0 ; proj < m_numProj ; proj++)
    {
        // loop on all X-ray energies
        for (size_t e = 0 ; e < m_energyVector.size() ; e++)
        {
            object->CreateObject(obj);
            diffuser->CreateDiffuser(diff);
            PropagateInFreeSpace();
            PropagateInMaterial();

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

            // must reset obj/diff matrices to zero to recreate them
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
