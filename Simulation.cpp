#include "Simulation.h"
#include "ImageExporter.h"
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
    m_M_obj.clear();
    m_M_diff.clear();

    for (int i = 0 ; i <= m_numObjVoxelsZ ; i++)
    {
        // distances are from the source to starting edge of the component
        // m_dM_obj.push_back(m_SDD/(m_SOD - m_objThickness/2 + (i - 1)*m_pixelSize/m_M_obj)); // dynamic magnification factors of slices
        // m_fM_obj.push_back((m_SOD - m_objThickness/2 + (i - 1)*m_pixelSize/m_M_obj)/(m_SOD - m_objThickness/2));
        m_M_obj.push_back(m_SDD/(m_SOD + i*m_pixelSize/m_SDD*m_SOD)); // dynamic magnification factors of slices
        m_fM_obj.push_back((m_SOD + i*m_pixelSize/m_SDD*m_SOD)/m_SOD); // Magnification after applying Fresnel scaling theory
    }

    for (int i = 0 ; i <= m_numDiffVoxelsZ ; i++)
    {
        // distances are from the source to starting edge of the component
        // m_dM_diff.push_back(m_SDD/(m_SdD - m_diffThickness/2 + (i - 1)*m_pixelSize/m_M_diff)); // dynamic magnification factors of slices
        // m_fM_diff.push_back((m_SdD - m_diffThickness/2 + (i - 1)*m_pixelSize/m_M_diff)/(m_SdD - m_diffThickness/2)); // Magnification after applying Fresnel scaling theory
        m_M_diff.push_back(m_SDD/(m_SdD + i*m_pixelSize/m_SDD*m_SdD)); // dynamic magnification factors of slices
        m_fM_diff.push_back((m_SdD + i*m_pixelSize/m_SDD*m_SdD)/m_SdD); // Magnification after applying Fresnel scaling theory
    }
}

void Simulation::setParams()
{    
    // tomo
    m_numProj = tab_setup->getNumProjections();

    // X-ray energy
    m_energyVector.clear(); // [keV]
    m_spectrumVector.clear();
    m_waveNumber.clear(); // [1/m]
    m_detEnergyResponse.clear();
    tab_source_detector->getXrayEnergy(m_energyVector, m_spectrumVector);
    m_waveNumber = tab_source_detector->getWaveNumber(m_energyVector);
    m_detEnergyResponse.assign(m_energyVector.size(), 1.0); // for now for simplicity

    // detector
    m_pixelSize = tab_source_detector->getPixelSizeMM();
    m_pixelSize = m_pixelSize*0.001; // [m]
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

    m_objThickness = tab_diffuser_object->getObjectThickness();
    m_diffThickness = tab_diffuser_object->getDiffuserThickness();

    m_numMVSlicesObj = tab_diffuser_object->getNumMVSlices(m_objThickness, m_numObjVoxelsZ, m_pixelSize, m_SDD/m_SOD);
    m_numMVSlicesDiff = tab_diffuser_object->getNumMVSlices(m_diffThickness, m_numDiffVoxelsZ, m_pixelSize, m_SDD/m_SdD);
    m_numInterp = tab_diffuser_object->getNumInterpolations();

    // mag factors based on the fresnel scaling theory
    getDynamicMagFactors();

    // // object type
    // if (tab_diffuser_object->m_cylinder){
    //     m_objectType = 1;

    // } else if (tab_diffuser_object->m_sphere){
    //     m_objectType = 2;

    // } else if (tab_diffuser_object->m_virtualObject){
    //     m_objectType = 3;
    // }

    // below part only for debugging
    // for (uint i = 0 ; i < m_energyVector.size() ; i++)
    // {
    //     qDebug() << "E : " << m_energyVector[i] << "\t" << m_spectrumVector[i];
    //     qDebug() << "k : " << m_waveNumber[i] << "\t" << m_waveNumber[i];
    // }

    // qDebug() << "pixel size : " << m_pixelSize;
    // qDebug() << "SOD : " << m_SOD;
    // qDebug() << "SdD : " << m_SdD;
    // qDebug() << "SDD : " << m_SDD;

    for (uint i = 0 ; i < m_physList.size() ; i++)
    {
        qDebug() << "physlist : d = " << m_physList[i].density << ", Z/A = " << m_physList[i].Z_A << ", name = " << m_physList[i].name;
    }

    qDebug() << "numPixels : " << m_numPixels;
    qDebug() << "obj-thickness : " << m_objThickness;
    qDebug() << "numVoxelsInZObj : " << m_numObjVoxelsZ;
    qDebug() << "numMVSlices_obj : " << m_numMVSlicesObj;
    qDebug() << "diff-thickness : " << m_diffThickness;
    qDebug() << "m_numMVSlicesDiff : " << m_numMVSlicesDiff;
    qDebug() << "m_numDiffVoxelsZ : " << m_numDiffVoxelsZ;
    qDebug() << "numInterp : " << m_numInterp;
    qDebug() << "numProj : " << m_numProj;
    qDebug() << "M_obj : " << m_M_obj;
    qDebug() << "M_diff : " << m_M_diff;
    qDebug() << "fM_obj : " << m_fM_obj;
    qDebug() << "fM_diff : " << m_fM_diff;
}

void Simulation::Gradient2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& Gx, std::vector<std::vector<double>>& Gy, const double spacing)
{
    int rows = matrix.size();
    int cols = matrix[0].size();

    // Compute the gradient in the x-direction (here between rows)
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 0; j < cols; ++j) {
            Gx[i][j] = (matrix[i + 1][j] - matrix[i - 1][j]) / (2.0 * spacing);
        }
    }
    // Handle boundary conditions for Gx (forward and backward difference)
    for (int j = 0; j < cols; ++j) {
        Gx[0][j] = (matrix[1][j] - matrix[0][j]) / spacing;  // Forward difference for top row
        Gx[rows - 1][j] = (matrix[rows - 1][j] - matrix[rows - 2][j]) / spacing;  // Backward difference for bottom row
    }

    // Compute the gradient in the y-direction (here between columns)
    for (int i = 0; i < rows; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            Gy[i][j] = (matrix[i][j + 1] - matrix[i][j - 1]) / (2.0 * spacing);
        }
    }
    // Handle boundary conditions for Gy (forward and backward difference)
    for (int i = 0; i < rows; ++i) {
        Gy[i][0] = (matrix[i][1] - matrix[i][0]) / spacing;  // Forward difference for left column
        Gy[i][cols - 1] = (matrix[i][cols - 1] - matrix[i][cols - 2]) / spacing;  // Backward difference for right column
    }

    // avoid nan values: not at this stage
    for (int i = 0 ; i < rows - 1 ; i++){
        for (int j = 0 ; j < cols -1 ; j++){
            if (std::isnan(Gx[i][j]) || std::isinf(Gx[i][j])){
                Gx[i][j] = 0;
            }
            if (std::isnan(Gy[i][j]) || std::isinf(Gy[i][j])){
                Gy[i][j] = 0;
            }
        }
    }

}

void Simulation::Laplacian2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& L, const double& spacing)
{
    int rows = matrix.size();
    int cols = matrix[0].size();

    // Compute the Laplacian for the interior points (central difference)
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            // Second derivative in the x-direction
            double d2x = (matrix[i + 1][j] - 2 * matrix[i][j] + matrix[i - 1][j]) / (spacing * spacing);

            // Second derivative in the y-direction
            double d2y = (matrix[i][j + 1] - 2 * matrix[i][j] + matrix[i][j - 1]) / (spacing * spacing);

            // Laplacian is the sum of second derivatives
            L[i][j] = d2x + d2y;
        }
    }

    // Handle boundary conditions using forward/backward differences
    // For the edges, we can use one-sided finite differences
    // Top and bottom rows
    for (int j = 0; j < cols; ++j) {
        L[0][j] = (matrix[1][j] - 2 * matrix[0][j] + matrix[0][j]) / (spacing * spacing); // forward difference
        L[rows - 1][j] = (matrix[rows - 1][j] - 2 * matrix[rows - 1][j] + matrix[rows - 2][j]) / (spacing * spacing); // backward difference
    }

    // Left and right columns
    for (int i = 0; i < rows; ++i) {
        L[i][0] = (matrix[i][1] - 2 * matrix[i][0] + matrix[i][0]) / (spacing * spacing); // forward difference
        L[i][cols - 1] = (matrix[i][cols - 1] - 2 * matrix[i][cols - 1] + matrix[i][cols - 2]) / (spacing * spacing); // backward difference
    }

    // avoid nan values
    for (int i = 0 ; i < rows - 1 ; i++){
        for (int j = 0 ; j < cols -1 ; j++){
            if (std::isnan(L[i][j]) || std::isinf(L[i][j])){
                L[i][j] = 0;
            }
        }
    }


}

 std::vector<std::vector<std::vector<int>>> Simulation::rotate3DObjectY(std::vector<std::vector<std::vector<int>>>& object, const int& numProjs)
{
    double rotAng = 360/numProjs * pi /180.0; // radians

    // rotation matrix around Y-axis
    double cosTheta = std::cos(rotAng);
    double sinTheta = std::sin(rotAng);

    // size of 3D object
    int Nx = object.size();
    int Ny = object[0].size();
    int Nz = object[0][0].size();

    std::vector<std::vector<std::vector<int>>> rotatedObject(Nx, std::vector<std::vector<int>>(Ny, std::vector<int>(Nz, 0)));

    // Get the center of the object (assume rotation around the center)
    int cx = Nx / 2;
    int cy = Ny / 2;
    int cz = Nz / 2;

    // Iterate through each point in the 3D object
    for (int x = 0; x < Nx; ++x) {
        for (int y = 0; y < Ny; ++y) {
            for (int z = 0; z < Nz; ++z) {

                // Translate point to the origin (center of object)
                int xt = x - cx;
                int yt = y - cy;
                int zt = z - cz;

                // Apply the rotation matrix (only Y-axis rotation here)
                double x_new = cosTheta * xt + sinTheta * zt;
                double y_new = yt;  // No change in y for Y-axis rotation
                double z_new = -sinTheta * xt + cosTheta * zt;

                // Translate back to the original coordinate system
                int xr = std::round(x_new + cx);
                int yr = std::round(y_new + cy);
                int zr = std::round(z_new + cz);

                // Check if the new coordinates are within bounds, if so, assign the value
                if (xr >= 0 && xr < Nx && yr >= 0 && yr < Ny && zr >= 0 && zr < Nz) {
                    rotatedObject[xr][yr][zr] = object[x][y][z];
                }
            }
        }
    }
    return rotatedObject;
}

double Simulation::Interpolation2D(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y,
                                   const std::vector<std::vector<std::vector<int>>>& grid, const int& sliceNum, const int& sliceStep, const double& n1_part,
                                   const double& n2_part, const double& mag, const double x, const double y, double default_value)
{
    int rows = grid.size();
    int cols = grid[0].size();

    // Find the four surrounding points in X and Y
    int x_low = 0, x_high = 0, y_low = 0, y_high = 0;

    // Find the correct x indices (x_low and x_high)
    for (int i = 0; i < rows - 1; i++) {
        if ( ((X[i][0]/mag) <= x) && ((X[i + 1][0]/mag) >= x) ) {
            x_low = i;
            x_high = i + 1;
            break;
        }
    }
    // Find the correct y indices (y_low and y_high)
    for (int i = 0; i < cols - 1; i++) {
        if ( ((Y[0][i]/mag) <= x) && ((Y[0][i + 1]/mag) >= x) ) {
            y_low = i;
            y_high = i + 1;
            break;
        }
    }

    // if ((grid[x_low][y_low][sliceNum+sliceStep] == 0) || (grid[x_low][y_low][sliceNum+sliceStep] == 0) || (grid[x_low][y_low][sliceNum+sliceStep] == 0) || (grid[x_low][y_low][sliceNum+sliceStep] == 0)){
    //     qDebug() << "x=" << x << ";y=" << y << ";x_low=" << x_low << ";x_high=" << x_high << ";y_low=" << y_low << ";y_high=" << y_high << "";
    //     qDebug() << "grid[x_low][y_low][sliceNum]=" << grid[x_low][y_low][sliceNum+sliceStep] << "grid[x_low][y_high][sliceNum]=" << grid[x_low][y_high][sliceNum+sliceStep];
    //     qDebug() << "grid[x_high][y_low][sliceNum]=" << grid[x_high][y_low][sliceNum+sliceStep] << "grid[x_high][y_high][sliceNum]=" << grid[x_high][y_high][sliceNum+sliceStep];
    // }
    // Check bounds, if out of range, return default value
    if (x_low < 0 || x_high >= rows || y_low < 0 || y_high >= cols || x_low == x_high || y_low == y_high) {
        return default_value;
        qDebug() << "default-value is used!";
    }

    // Bilinear interpolation
    double x_ratio = (x - X[x_low][0]/mag) / (X[x_high][0]/mag - X[x_low][0]/mag);
    double y_ratio = (y - Y[0][y_low]/mag) / (Y[0][y_high]/mag - Y[0][y_low]/mag);

    double top = ( (grid[x_low][y_low][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_low][y_low][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_low][y_low][sliceNum+sliceStep]) ) * (1 - y_ratio)
                 + ( (grid[x_low][y_high][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_low][y_high][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_low][y_high][sliceNum+sliceStep]) ) * y_ratio;

    double bottom = ( (grid[x_high][y_low][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_high][y_low][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_high][y_low][sliceNum+sliceStep]) ) * (1 - y_ratio)
                    + ( (grid[x_high][y_high][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_high][y_high][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_high][y_high][sliceNum+sliceStep]) ) * y_ratio;

        return ((1 - x_ratio) * top + x_ratio * bottom); // still doesn't handle point 0,0, even with numInterp=1 : TODO

}

void Simulation::SubPixelInterpolation(const std::vector<std::vector<std::vector<int>>>& subject, std::vector<std::vector<double>>& beta,
                                       std::vector<std::vector<double>>& delta, const Materials::refractiveIndex& n1, const Materials::refractiveIndex& n2,
                                       const double& sourceToSubjectDist, const std::vector<double>& Mag, const int& sliceNum, const int& sliceStep, const int& m_numInterp)
{
    double M_interp = 0;
    beta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // imaginary part of the complex refractive index
    delta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // real part of the complex refractive index

    // loop over different interpolation steps and accumulating them
    // for (int interp = 0 ; interp < m_numInterp ; interp++)
    // {
        // M_interp = m_SDD/(sourceToSubjectDist - thickness/2 + (sliceNum-1)*m_pixelSize/globalMag + interp*m_pixelSize/globalMag/m_numInterp);
        // M_interp = m_SDD/(sourceToSubjectDist + sliceNum*m_pixelSize/globalMag + interp*m_pixelSize/globalMag/m_numInterp);

        // M_interp = m_SDD/(sourceToSubjectDist + sliceNum*m_pixelSize/Mag[0] + (interp/m_numInterp)*(sliceStep)*m_pixelSize/Mag[0]);
        // qDebug() << "M_interp=" <<  M_interp;

        for (int i = 0 ; i < m_numPixels ; i++){
            for (int j = 0 ; j < m_numPixels ; j++){
                beta[i][j] = ( (subject[i][j][sliceNum] == 1) ? ((n1.imaginaryPart) * subject[i][j][sliceNum]) : ((n2.imaginaryPart/2) * subject[i][j][sliceNum]) );
                delta[i][j] = ( (subject[i][j][sliceNum] == 1) ? ((n1.realPart) * subject[i][j][sliceNum]) : ((n2.realPart/2) * subject[i][j][sliceNum]) );
                // beta[i][j]+=Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.imaginaryPart, n2.imaginaryPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
                // delta[i][j]+=Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.realPart, n2.realPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
                // if ((i==5 && j==5) || (i==255 && j==2)){
                    // qDebug() << "interp #" << (interp + 1) << " : Interp2D_beta(" << i << "," << j << ")=" << Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.imaginaryPart, n2.imaginaryPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
                    // qDebug() << "interp #" << (interp + 1) << " : Interp2D_delta(" << i << "," << j << ")=" << Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.realPart, n2.realPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
                // }
            }
        }
        // qDebug() << "interpolation #" << (interp + 1);
    // }

    // // averaging over interpolations
    // for (int i = 0 ; i < m_numPixels ; i++){
    //     for (int j = 0 ; j < m_numPixels ; j++){
    //         beta[i][j] = beta[i][j] / m_numInterp;
    //         delta[i][j]= delta[i][j] / m_numInterp;
    //     }
    // }
}

void Simulation::PropagateInMaterial(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const std::vector<std::vector<std::vector<int>>>& subject, const double& thickness,
                                     const double& sourceToSubjectDist, const std::vector<double>& Mag, const std::vector<double>& fresnelMag,const int& numMVSlices,
                                     const int& numVoxelSlicesInZ, const Materials::refractiveIndex& refrIndx1, const Materials::refractiveIndex& refrIndx2, const size_t& energyIndx)
{
    // declarations for each slice
    double coneBeamSliceThickness = 0; // cone-beam slice thickness at plane
    int sliceStep = numVoxelSlicesInZ/numMVSlices;
    qDebug() << "refrIndx1=" << &refrIndx1 << "; refrIndx2=" << &refrIndx2;

    // create magnified thickness mesh and interpolate refrative index in sub-pixel resolution
    for (int s = 0 ; s < numVoxelSlicesInZ ; s+=sliceStep)
    {
        Gx_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
        Gy_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
        Gx_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
        Gy_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
        L_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));

        // double GxIGxphi = 0, GyIGyphi = 0, thirdTerm = 0;

        SubPixelInterpolation(subject, beta, delta, refrIndx1, refrIndx2, sourceToSubjectDist, Mag, s, sliceStep, m_numInterp); // sub-pixel interpolation of the complex refrative index
        Gradient2D(I, Gx_I, Gy_I, m_pixelSize/Mag[s+sliceStep]);
        Gradient2D(phi, Gx_phi, Gy_phi, m_pixelSize/Mag[s+sliceStep]);
        Laplacian2D(phi, L_phi, m_pixelSize/Mag[s+sliceStep]);

        // qDebug() << "m_waveNumber[energyIndx]" << m_waveNumber[energyIndx];
        // qDebug() << "sliceNum=" << s+sliceStep;


        // update cone-beam magnified slice thickness and then I and phi
        for (int i = 0 ; i < m_numPixels ; i++){
            for (int j = 0 ; j < m_numPixels ; j++){

                // coneBeamSliceThickness = std::sqrt( m_rsqr[i][j]/(std::pow(dynamicMag[s], 2)) + std::pow((sourceToSubjectDist - thickness/2 + s*m_pixelSize/globalMag), 2) )
                //                         - std::sqrt( m_rsqr[i][j]/(std::pow(dynamicMag[s-sliceStep-1], 2)) + std::pow((sourceToSubjectDist - thickness/2 + (s-sliceStep-1)*m_pixelSize/globalMag), 2) );
                // coneBeamSliceThickness = std::sqrt( m_rsqr[i][j]/(std::pow(dynamicMag[s+sliceStep], 2)) + std::pow((sourceToSubjectDist + (s+sliceStep)*m_pixelSize/globalMag), 2) )
                //                          - std::sqrt( m_rsqr[i][j]/(std::pow(dynamicMag[s], 2)) + std::pow((sourceToSubjectDist + (s)*m_pixelSize/globalMag), 2) );
                coneBeamSliceThickness = std::sqrt( m_rsqr[i][j]/(std::pow(Mag[s+sliceStep], 2)) + std::pow((sourceToSubjectDist + (s+sliceStep)*m_pixelSize/Mag[0]), 2) )
                                         - std::sqrt( m_rsqr[i][j]/(std::pow(Mag[s], 2)) + std::pow((sourceToSubjectDist + (s)*m_pixelSize/Mag[0]), 2) );
                // loop to update I and phi

                // if ( (i==5 && j==5) || (i==255 && j==2)){

                // if (std::isinf((Gx_I[i][j] * Gx_phi[i][j])) || std::isnan((Gx_I[i][j] * Gx_phi[i][j]))) {
                        // qDebug() << "Gx_I[i][j] * Gx_phi[i][j] is inf or nan in (" << i << "," << j << ")";
                        // qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
                        // qDebug() << "coneBeamSliceThickness = " << coneBeamSliceThickness;
                        // qDebug() << "(i, j)=(" << i << "," << j << ") -> beta[i][j]=" << beta[i][j] << " ; delta[i][j]=" << delta[i][j];
                        // qDebug() << "GxI = " << Gx_I[i][j] << ", GyI = " << Gy_I[i][j] << ", Gxphi = " << Gx_phi[i][j] << ", Gyphi = " << Gy_phi[i][j] << ", Lphi = " << L_phi[j][j];
                        // qDebug() << "(I-a*b)*c, a->" << (m_pixelSize/Mag[0]/fresnelMag[s+sliceStep]/m_waveNumber[energyIndx]);
                        // qDebug() << "(I-a*b)*c, b->" << (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]);
                        // qDebug() << "(I-a*b)*c, c->" << std::exp(-2 * m_waveNumber[energyIndx] * beta[i][j] * coneBeamSliceThickness);
                        // qDebug() << "phi-d, d->" << m_waveNumber[energyIndx] * delta[i][j] * coneBeamSliceThickness;
                        // qDebug() << "I=" << I[i][j];
                        // qDebug() << "phi=" << phi[i][j];
                        // qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
                //     } else if (std::isinf((Gy_I[i][j] * Gy_phi[i][j])) || std::isnan((Gy_I[i][j] * Gy_phi[i][j]))) {
                //         qDebug() << "Gy_I[i][j] * Gy_phi[i][j] is inf or nan in (" << i << "," << j << ")";
                //         qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
                //         qDebug() << "coneBeamSliceThickness = " << coneBeamSliceThickness;
                //         qDebug() << "(i, j)=(" << i << "," << j << ") -> beta[i][j]=" << beta[i][j] << " ; delta[i][j]=" << delta[i][j];
                //         qDebug() << "GxI = " << Gx_I[i][j] << ", GyI = " << Gy_I[i][j] << ", Gxphi = " << Gx_phi[i][j] << ", Gyphi = " << Gy_phi[i][j] << ", Lphi = " << L_phi[j][j];
                //         qDebug() << "(I-a*b)*c, a->" << (m_pixelSize/Mag[0]/fresnelMag[s+sliceStep]/m_waveNumber[energyIndx]);
                //         qDebug() << "(I-a*b)*c, b->" << (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]);
                //         qDebug() << "(I-a*b)*c, c->" << std::exp(-2 * m_waveNumber[energyIndx] * beta[i][j] * coneBeamSliceThickness);
                //         qDebug() << "phi-d, d->" << m_waveNumber[energyIndx] * delta[i][j] * coneBeamSliceThickness;
                //         qDebug() << "I=" << I[i][j];
                //         qDebug() << "phi=" << phi[i][j];
                //         qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
                //     } else if (std::isinf((4 * I[i][j] * L_phi[i][j])) || std::isnan((4 * I[i][j] * L_phi[i][j]))) {
                //         qDebug() << "4 * I[i][j] * L_phi[i][j] is inf or nan in (" << i << "," << j << ")";
                //         qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
                //         qDebug() << "coneBeamSliceThickness = " << coneBeamSliceThickness;
                //         qDebug() << "(i, j)=(" << i << "," << j << ") -> beta[i][j]=" << beta[i][j] << " ; delta[i][j]=" << delta[i][j];
                //         qDebug() << "GxI = " << Gx_I[i][j] << ", GyI = " << Gy_I[i][j] << ", Gxphi = " << Gx_phi[i][j] << ", Gyphi = " << Gy_phi[i][j] << ", Lphi = " << L_phi[j][j];
                //         qDebug() << "(I-a*b)*c, a->" << (m_pixelSize/Mag[0]/fresnelMag[s+sliceStep]/m_waveNumber[energyIndx]);
                //         qDebug() << "(I-a*b)*c, b->" << (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]);
                //         qDebug() << "(I-a*b)*c, c->" << std::exp(-2 * m_waveNumber[energyIndx] * beta[i][j] * coneBeamSliceThickness);
                //         qDebug() << "phi-d, d->" << m_waveNumber[energyIndx] * delta[i][j] * coneBeamSliceThickness;
                //         qDebug() << "I=" << I[i][j];
                //         qDebug() << "phi=" << phi[i][j];
                //         qDebug() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
                //     }

                // }

                I[i][j] = ( I[i][j] - (m_pixelSize/Mag[0]/fresnelMag[s+sliceStep]/m_waveNumber[energyIndx]) *
                                         (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) ) *
                            (std::exp(-2 * m_waveNumber[energyIndx] * beta[i][j] * coneBeamSliceThickness));

                phi[i][j] = phi[i][j] - (m_waveNumber[energyIndx] * delta[i][j] * coneBeamSliceThickness);

                // if ( (i==5 && j==5) || (i==255 && j==2) ){

                    // qDebug() << "I=" << I[i][j];
                    // qDebug() << "phi=" << phi[i][j];
                    // qDebug() << "========================================================================================";
                // }
            }
        }
        // qDebug() << "========================================================================================";
    }
}

void Simulation::PropagateInFreeSpace(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const double& thickness1, const double& dist1, const double& dist2,
                                      const double& mag2, const double& mag1, const size_t& energyIndx, const bool isLastPropagation)
{
    // 1: first located element
    // thickness1 or 2: thickness of the first or second component from the source : i.e., object, diffuser, or detector(in this case=0)
    // dist1 or 2 : source to 1st or 2nd component distance : i.e. SOD, SDD, or SDD
    // Mag1 or 2 : relevant magnification factor(global/dynamic) of the 1st or 2nd located component (from the source)

    double energyWeight, detectorResponse = 0;
    double propagDist = 0;
    Gx_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
    Gy_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
    Gx_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
    Gy_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
    L_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));

    Gradient2D(I, Gx_I, Gy_I, m_pixelSize/mag1);
    Gradient2D(phi, Gx_phi, Gy_phi, m_pixelSize/mag1);
    Laplacian2D(phi, L_phi, m_pixelSize/mag1);

    if (isLastPropagation){
        energyWeight = m_spectrumVector[energyIndx];
        detectorResponse = m_detEnergyResponse[energyIndx];
        qDebug() << "It is the last propagation";
    } else {
        energyWeight = 1;
        detectorResponse = 1;
    }

    for (int i = 0 ; i < m_numPixels ; i++){
        for (int j = 0 ; j< m_numPixels ; j++){

            // propagDist = std::sqrt(std::pow((dist2 - thickness2/2), 2) + m_rsqr[i][j]/(globalMag2*globalMag2))
            //              - std::sqrt(std::pow((dist1+thickness1/2),2) + m_rsqr[i][j]/(dynamicMag1*dynamicMag1));
            propagDist = std::sqrt(std::pow(dist2, 2) + m_rsqr[i][j]/(mag2*mag2))
                         - std::sqrt(std::pow((dist1+thickness1),2) + m_rsqr[i][j]/(mag1*mag1));

            I[i][j] = energyWeight * detectorResponse * ( I[i][j] - (propagDist/dist2*dist1/m_waveNumber[energyIndx]) * (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) );

            // if ( (i==5 && j==5) || (i==255 && j==2) ){
            //     qDebug() << "****************************************";
            //     qDebug() << "propagdist=" << propagDist;
            //     qDebug() << "I=" << I[i][j];
            //     qDebug() << "phi=" << phi[i][j];
            //     qDebug() << "****************************************";
            // }
        }
    }
}

void Simulation::ConeBeamSimulation()
{
    // temporary matrices in next loops
    beta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // imaginary part of the complex refractive index
    delta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // real part of the complex refractive index
    Gx_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // intensity gradient along x : here vertical side
    Gy_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // intensity gradient along y : here horizontal side
    Gx_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // phase gradient along x : here vertical side
    Gy_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // phase gradient along y : here horizontal side
    L_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // 2D Laplacian of phi

    setParams(); // set all parameters

    // object declarations
    std::vector<std::vector<std::vector<int>>> obj(m_numPixels, std::vector<std::vector<int>>(m_numPixels, std::vector<int>(m_numObjVoxelsZ, 0)));
    std::vector<Materials::refractiveIndex> n_obj = mat->RefractiveIndex(m_energyVector, m_physList[3].formula, m_physList[3].density);
    std::unique_ptr<Object> object = std::make_unique<Object>(m_numObjVoxelsZ, m_numPixels, m_pixelSize/m_M_obj[0], tab_diffuser_object);

    // diffuser declarations
    std::vector<std::vector<std::vector<int>>> diff(m_numPixels, std::vector<std::vector<int>>(m_numPixels, std::vector<int>(m_numDiffVoxelsZ, 2)));
    std::vector<Materials::refractiveIndex> n_diff = mat->RefractiveIndex(m_energyVector, m_physList[1].formula, m_physList[1].density);
    std::vector<Materials::refractiveIndex> n_base = mat->RefractiveIndex(m_energyVector, m_physList[2].formula, m_physList[2].density);
    int numGrits = 0;
    std::unique_ptr<Diffuser> diffuser = std::make_unique<Diffuser>(m_numPixels, m_numDiffVoxelsZ, m_pixelSize/m_M_diff[0], tab_diffuser_object, numGrits);
    qDebug() << "numGrits=" << numGrits;

    // loop over projections
    for (int proj = 0 ; proj < m_numProj ; proj++)
    {
        // set intensity distribution's amplitude and phase (before object or diffuser)
        m_I_bg.assign(m_numPixels, std::vector<double>(m_numPixels, 1.0));
        m_phi_bg.assign(m_numPixels, std::vector<double>(m_numPixels, 0.0));
        m_I_fg.assign(m_numPixels, std::vector<double>(m_numPixels, 1.0));
        m_phi_fg.assign(m_numPixels, std::vector<double>(m_numPixels, 0.0));

        // loop on all X-ray energies
        for (size_t e = 0 ; e < m_energyVector.size() ; e++)
        {
            object->CreateObject(obj); // create object
            diffuser->CreateDiffuser(diff); // create diffuser

            // qDebug() << "diff[120][130]=" << diff[120][130] << "; diff[2][2]=" << diff[2][2] << "; diff[120][2]=" << diff[120][2] << "; diff[2][130]=" << diff[2][130] << "\n";
            // qDebug() << "obj[120][130]=" << obj[120][130] << "; obj[2][2]=" << obj[2][2] << "; obj[3][3]=" << obj[3][3] << "\n";
            // qDebug() << "diff[5][5]=" << diff[5][5] << "; diff[255][2]=" << diff[255][2] << "\n";
            // qDebug() << "obj[5][5]=" << obj[5][5] << "; obj[255][2]=" << obj[255][2] << "\n";

            qDebug() << "\nn_obj-im : " << n_obj[e].imaginaryPart << "; n_obj-re : " << n_obj[e].realPart;
            qDebug() << "n_diff-im : " << n_diff[e].imaginaryPart << "; n_diff-re : " << n_diff[e].realPart;
            qDebug() << "n_base-im : " << n_base[e].imaginaryPart << "; n_base-re : " << n_base[e].realPart;

            // reference image (only with diffuser)
            qDebug() << "only diff-> propagate in diff :\n";
            PropagateInMaterial(m_I_bg, m_phi_bg, diff, m_diffThickness, m_SdD, m_M_diff, m_fM_diff, m_numMVSlicesDiff, m_numDiffVoxelsZ, n_diff[e], n_base[e], e); // propagate through the diffuser
            // qDebug() << "\nm_I_bg[5][5]=" << m_I_bg[5][5];
            // qDebug() << "m_phi_bg[5][5]=" << m_phi_bg[5][5];
            // qDebug() << "m_I_bg[255][2]=" << m_I_bg[255][2];
            // qDebug() << "m_phi_bg[255][2]=" << m_phi_bg[255][2];
            qDebug() << "\npropagate from diff to det:\n";
            PropagateInFreeSpace(m_I_bg, m_phi_bg, m_diffThickness, m_SdD, m_SDD, 1, m_M_diff[m_numDiffVoxelsZ], e, true); // propagate from diffuser to detector
            // qDebug() << "\nm_I_bg[5][5]=" << m_I_bg[5][5];
            // qDebug() << "m_phi_bg[5][5]=" << m_phi_bg[5][5];
            // qDebug() << "m_I_bg[255][2]=" << m_I_bg[255][2];
            // qDebug() << "m_phi_bg[255][2]=" << m_phi_bg[255][2];


            // image (diffuser + object)
            if (m_SOD < m_SdD){ // object is before diffuser

                qDebug() << "with diff&obj->propagate in object:\n";
                PropagateInMaterial(m_I_fg, m_phi_fg, obj, m_objThickness, m_SOD, m_M_obj, m_fM_obj, m_numMVSlicesObj, m_numObjVoxelsZ, n_obj[e], n_obj[e], e); // propagate through the object
                // qDebug() << "\nm_I_fg[5][5]=" << m_I_fg[5][5];
                // qDebug() << "m_phi_fg[5][5]=" << m_phi_fg[5][5];
                // qDebug() << "m_I_fg[255][2]=" << m_I_fg[255][2];
                // qDebug() << "m_phi_fg[255][2]=" << m_phi_fg[255][2];
                qDebug() << "\npropagate from obj to diff:\n";
                PropagateInFreeSpace(m_I_fg, m_phi_fg, m_objThickness, m_SOD, m_SdD, m_M_diff[0], m_M_obj[m_numObjVoxelsZ], e, false); // propagate from object to diffuser (or from diffuser to object)
                // qDebug() << "\nm_I_fg[5][5]=" << m_I_fg[5][5];
                // qDebug() << "m_phi_fg[5][5]=" << m_phi_fg[5][5];
                // qDebug() << "m_I_fg[255][2]=" << m_I_fg[255][2];
                // qDebug() << "m_phi_fg[255][2]=" << m_phi_fg[255][2];
                qDebug() << "\npropagate in diff:\n";
                PropagateInMaterial(m_I_fg, m_phi_fg, diff, m_diffThickness, m_SdD, m_M_diff, m_fM_diff,m_numMVSlicesDiff, m_numDiffVoxelsZ, n_diff[e], n_base[e], e); // propagate through the diffuser
                // qDebug() << "\nm_I_fg[5][5]=" << m_I_fg[5][5];
                // qDebug() << "m_phi_fg[5][5]=" << m_phi_fg[5][5];
                // qDebug() << "m_I_fg[255][2]=" << m_I_fg[255][2];
                // qDebug() << "m_phi_fg[255][2]=" << m_phi_fg[255][2];
                qDebug() << "\npropagate from diff to det:\n";
                PropagateInFreeSpace(m_I_fg, m_phi_fg, m_diffThickness, m_SdD, m_SDD, 1, m_M_diff[m_numDiffVoxelsZ], e, true); // propagate from diffuser to detector
                // qDebug() << "\nm_I_fg[5][5]=" << m_I_fg[5][5];
                // qDebug() << "m_phi_fg[5][5]=" << m_phi_fg[5][5];
                // qDebug() << "m_I_fg[255][2]=" << m_I_fg[255][2];
                // qDebug() << "m_phi_fg[255][2]=" << m_phi_fg[255][2];

            } else { // object is after diffuser

                qDebug() << "obj+diff : propagate in diff :";
                PropagateInMaterial(m_I_fg, m_phi_fg, diff, m_diffThickness, m_SdD, m_M_diff, m_fM_diff, m_numMVSlicesDiff, m_numDiffVoxelsZ, n_diff[e], n_base[e], e); // propagate through the diffuser
                qDebug() << "propagate from diff to obj:";
                PropagateInFreeSpace(m_I_fg, m_phi_fg, m_diffThickness, m_SdD, m_SOD, m_M_obj[0], m_M_diff[m_numDiffVoxelsZ], e, false); // propagate from diffuser to object
                qDebug() << "propagate in obj :";
                PropagateInMaterial(m_I_fg, m_phi_fg, obj, m_objThickness, m_SOD, m_M_obj, m_fM_obj, m_numMVSlicesObj, m_numObjVoxelsZ, n_obj[e], n_obj[e], e); // propagate through the object
                qDebug() << "propagate from obj to det:";
                PropagateInFreeSpace(m_I_fg, m_phi_fg, m_objThickness, m_SOD, m_SDD, 1, m_M_obj[m_numObjVoxelsZ], e, true); // propagate from object to detector

            }


        }

        // qDebug() << "\nm_I_bg[5][5]=" << m_I_bg[5][5];
        // qDebug() << "m_phi_bg[5][5]=" << m_phi_bg[5][5];
        // qDebug() << "m_I_bg[255][2]=" << m_I_bg[255][2];
        // qDebug() << "m_phi_bg[255][2]=" << m_phi_bg[255][2];
        // qDebug() << "\nm_I_fg[5][5]=" << m_I_fg[5][5];
        // qDebug() << "m_phi_fg[5][5]=" << m_phi_fg[5][5];
        // qDebug() << "m_I_fg[255][2]=" << m_I_fg[255][2];
        // qDebug() << "m_phi_fg[255][2]=" << m_phi_fg[255][2];

        // int c1=0, c2=0, c3=0, c4 = 0;

        // for (size_t  i = 0 ; i < m_I_bg.size() ; i++){
        //     for (size_t  j = 0 ; j < m_I_bg[0].size() ; j++){
        //         if (std::isnan(m_I_bg[i][j]) || std::isinf(m_I_bg[i][j])){
        //             qDebug() << "m_I_bg(" << i << "," << j << ")=" << m_I_bg[i][j];
        //             c1++;
        //         } else if (std::isnan(m_phi_bg[i][j]) || std::isinf(m_phi_bg[i][j])){
        //             qDebug() << "m_phi_bg(" << i << "," << j << ")=" << m_phi_bg[i][j];
        //             c2++;
        //         } else if (std::isnan(m_I_fg[i][j]) || std::isinf(m_I_fg[i][j])){
        //             qDebug() << "m_I_fg(" << i << "," << j << ")=" << m_I_fg[i][j];
        //             c3++;
        //         } else if (std::isnan(m_phi_fg[i][j]) || std::isinf(m_phi_fg[i][j])){
        //             qDebug() << "m_phi_fg(" << i << "," << j << ")=" << m_phi_fg[i][j];
        //             c4++;
        //         }
        //     }
        // }
        // qDebug() << "c1=" << c1 << ", c2=" << c2 << ", c3=" << c3 << ", c4=" << c4 <<" pixels are equal nan or inf.";

        // Apply detector PSF to final images
        // BlurImage();

        // TODO : user determines the output image name and format
        // save background (bg) and foreground (fg) images to files
        qDebug() << "writing images to files ... :";
        // std::string output_image_bg = "Image_bg_" + std::to_string(proj+1) + ".nii";
        // std::string output_image_fg = "Image_fg_" + std::to_string(proj+1) + ".nii";

        // std::unique_ptr<ImageExporter<double>> image_bg = std::make_unique<ImageExporter<double>>(static_cast<size_t>(m_numPixels), static_cast<size_t>(m_numPixels));
        // image_bg->SetData(m_I_bg);
        // image_bg->Save2D(output_image_bg);

        // std::unique_ptr<ImageExporter<double>> image_fg = std::make_unique<ImageExporter<double>>(static_cast<size_t>(m_numPixels), static_cast<size_t>(m_numPixels));
        // image_fg->SetData(m_I_fg);
        // image_fg->Save2D(output_image_fg);

        // // rotate object
        // obj = rotate3DObjectY(obj, m_numProj);


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

    }
}

void Simulation::ParBeamSimulation()
{

}

void Simulation::InitiateSimulation()
{
    if (tab_source_detector->m_parBeam){
        ParBeamSimulation();

    } else if (tab_source_detector->m_coneBeam){
        ConeBeamSimulation();

    } else {
        QMessageBox::warning(this, "ERROR", "Please select either parallel beam or cone-beam geometry. Aborting ...");
        return;

    }
}

//TODO
// void Simulation::WriteToImage()
// {

// }
