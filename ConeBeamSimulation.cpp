// #include "ConeBeamSimulation.h"
// #include "Object.h"
// #include "Diffuser.h"

// #include <QMessageBox>
// #include <QDebug>



// ConeBeamSimulation::ConeBeamSimulation(Simulation* sim) : simPtr(sim)//, compute(new Compute())
// {
//     // compute = new Compute();

//     m_numPixels = simPtr->m_numPixels;
//     m_numObjVoxelsZ = simPtr->m_numObjVoxelsZ;
//     m_numMVSlicesObj = simPtr->m_numMVSlicesObj;
//     m_numMVSlicesDiff = simPtr->m_numMVSlicesDiff;
//     m_numDiffVoxelsZ = simPtr->m_numDiffVoxelsZ;
//     m_numInterp = simPtr->m_numInterp;
//     m_numProj = simPtr->m_numProj;
//     m_pixelSize = simPtr->m_pixelSize;
//     m_objThickness = simPtr->m_objThickness;
//     m_diffThickness = simPtr->m_diffThickness;
//     m_SOD = simPtr->m_SDD;
//     m_SdD = simPtr->m_SdD;
//     m_SDD = simPtr->m_SDD;

// }

// ConeBeamSimulation::~ConeBeamSimulation()
// {
//     // delete compute;
//     delete simPtr;
// }


// // void ConeBeamSimulation::Gradient2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& Gx, std::vector<std::vector<double>>& Gy, const double spacing)
// // {
// //     int rows = matrix.size();
// //     int cols = matrix[0].size();

// //     // Compute the gradient in the x-direction (here between rows)
// //     for (int i = 1; i < rows - 1; ++i) {
// //         for (int j = 0; j < cols; ++j) {
// //             Gx[i][j] = (matrix[i + 1][j] - matrix[i - 1][j]) / (2.0 * spacing);
// //         }
// //     }
// //     // Handle boundary conditions for Gx (forward and backward difference)
// //     for (int j = 0; j < cols; ++j) {
// //         Gx[0][j] = (matrix[1][j] - matrix[0][j]) / spacing;  // Forward difference for top row
// //         Gx[rows - 1][j] = (matrix[rows - 1][j] - matrix[rows - 2][j]) / spacing;  // Backward difference for bottom row
// //     }

// //     // Compute the gradient in the y-direction (here between columns)
// //     for (int i = 0; i < rows; ++i) {
// //         for (int j = 1; j < cols - 1; ++j) {
// //             Gy[i][j] = (matrix[i][j + 1] - matrix[i][j - 1]) / (2.0 * spacing);
// //         }
// //     }
// //     // Handle boundary conditions for Gy (forward and backward difference)
// //     for (int i = 0; i < rows; ++i) {
// //         Gy[i][0] = (matrix[i][1] - matrix[i][0]) / spacing;  // Forward difference for left column
// //         Gy[i][cols - 1] = (matrix[i][cols - 1] - matrix[i][cols - 2]) / spacing;  // Backward difference for right column
// //     }

// //     // avoid nan values: not at this stage
// //     for (int i = 0 ; i < rows - 1 ; i++){
// //         for (int j = 0 ; j < cols -1 ; j++){
// //             if (std::isnan(Gx[i][j]) || std::isinf(Gx[i][j])){
// //                 Gx[i][j] = 0;
// //             }
// //             if (std::isnan(Gy[i][j]) || std::isinf(Gy[i][j])){
// //                 Gy[i][j] = 0;
// //             }
// //         }
// //     }

// // }

// // void ConeBeamSimulation::Laplacian2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& L, const double& spacing)
// // {
// //     int rows = matrix.size();
// //     int cols = matrix[0].size();

// //     // Compute the Laplacian for the interior points (central difference)
// //     for (int i = 1; i < rows - 1; ++i) {
// //         for (int j = 1; j < cols - 1; ++j) {
// //             // Second derivative in the x-direction
// //             double d2x = (matrix[i + 1][j] - 2 * matrix[i][j] + matrix[i - 1][j]) / (spacing * spacing);

// //             // Second derivative in the y-direction
// //             double d2y = (matrix[i][j + 1] - 2 * matrix[i][j] + matrix[i][j - 1]) / (spacing * spacing);

// //             // Laplacian is the sum of second derivatives
// //             L[i][j] = d2x + d2y;
// //         }
// //     }

// //     // Handle boundary conditions using forward/backward differences
// //     // For the edges, we can use one-sided finite differences
// //     // Top and bottom rows
// //     for (int j = 0; j < cols; ++j) {
// //         L[0][j] = (matrix[1][j] - 2 * matrix[0][j] + matrix[0][j]) / (spacing * spacing); // forward difference
// //         L[rows - 1][j] = (matrix[rows - 1][j] - 2 * matrix[rows - 1][j] + matrix[rows - 2][j]) / (spacing * spacing); // backward difference
// //     }

// //     // Left and right columns
// //     for (int i = 0; i < rows; ++i) {
// //         L[i][0] = (matrix[i][1] - 2 * matrix[i][0] + matrix[i][0]) / (spacing * spacing); // forward difference
// //         L[i][cols - 1] = (matrix[i][cols - 1] - 2 * matrix[i][cols - 1] + matrix[i][cols - 2]) / (spacing * spacing); // backward difference
// //     }

// //     // avoid nan values
// //     for (int i = 0 ; i < rows - 1 ; i++){
// //         for (int j = 0 ; j < cols -1 ; j++){
// //             if (std::isnan(L[i][j]) || std::isinf(L[i][j])){
// //                 L[i][j] = 0;
// //             }
// //         }
// //     }


// // }

// // std::vector<std::vector<std::vector<int>>> ConeBeamSimulation::rotate3DObjectY(std::vector<std::vector<std::vector<int>>>& object, const int& numProjs)
// // {
// //     double rotAng = 360/numProjs * pi /180.0; // radians

// //     // rotation matrix around Y-axis
// //     double cosTheta = std::cos(rotAng);
// //     double sinTheta = std::sin(rotAng);

// //     // size of 3D object
// //     int Nx = object.size();
// //     int Ny = object[0].size();
// //     int Nz = object[0][0].size();

// //     std::vector<std::vector<std::vector<int>>> rotatedObject(Nx, std::vector<std::vector<int>>(Ny, std::vector<int>(Nz, 0)));

// //     // Get the center of the object (assume rotation around the center)
// //     int cx = Nx / 2;
// //     int cy = Ny / 2;
// //     int cz = Nz / 2;

// //     // Iterate through each point in the 3D rotated object -> backward mapping
// //     for (int xr = 0; xr < Nx; ++xr) {
// //         for (int yr = 0; yr < Ny; ++yr) {
// //             for (int zr = 0; zr < Nz; ++zr) {

// //                 // Translate point to the origin (center of object)
// //                 int xt = xr - cx;
// //                 int yt = yr - cy;
// //                 int zt = zr - cz;

// //                 // Apply the rotation matrix (only Y-axis rotation here)
// //                 double xo = cosTheta * xt - sinTheta * zt;
// //                 double yo = yt;  // No change in y for Y-axis rotation
// //                 double zo = sinTheta * xt + cosTheta * zt;

// //                 // Translate back to the original coordinate system
// //                 int xo_int = std::round(xo + cx);
// //                 int yo_int = std::round(yo + cy);
// //                 int zo_int = std::round(zo + cz);

// //                 // Check if the new coordinates are within bounds, if so, assign the value
// //                 if (xo_int >= 0 && xo_int < Nx && yo_int >= 0 && yo_int < Ny && zo_int >= 0 && zo_int < Nz) {
// //                     rotatedObject[xr][yr][zr] = object[xo_int][yo_int][zo_int];
// //                 }
// //             }
// //         }
// //     }
// //     return rotatedObject;
// // }

// // double ConeBeamSimulation::Interpolation2D(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y,
// //                                 const std::vector<std::vector<std::vector<int>>>& grid, const int& sliceNum, const int& sliceStep, const double& n1_part,
// //                                 const double& n2_part, const double& mag, const double x, const double y, double default_value)
// // {
// //     int rows = grid.size();
// //     int cols = grid[0].size();

// //     double x_ratio = 0, y_ratio = 0;

// //     // Find the four surrounding points in X and Y
// //     int x_low = 0, x_high = 0, y_low = 0, y_high = 0;

// //     // Find the correct x indices (x_low and x_high)
// //     for (int i = 0; i < rows - 1; i++) {
// //         if ( ((X[i][0]/mag) <= x) && ((X[i + 1][0]/mag) >= x) ) {
// //             x_low = i;
// //             x_high = i + 1;
// //             break;
// //         }
// //     }
// //     // Find the correct y indices (y_low and y_high)
// //     for (int i = 0; i < cols - 1; i++) {
// //         if ( ((Y[0][i]/mag) <= x) && ((Y[0][i + 1]/mag) >= x) ) {
// //             y_low = i;
// //             y_high = i + 1;
// //             break;
// //         }
// //     }

// //     // if ((grid[x_low][y_low][sliceNum+sliceStep] == 0) || (grid[x_low][y_low][sliceNum+sliceStep] == 0) || (grid[x_low][y_low][sliceNum+sliceStep] == 0) || (grid[x_low][y_low][sliceNum+sliceStep] == 0)){
// //     //     qDebug() << "x=" << x << ";y=" << y << ";x_low=" << x_low << ";x_high=" << x_high << ";y_low=" << y_low << ";y_high=" << y_high << "";
// //     //     qDebug() << "grid[x_low][y_low][sliceNum]=" << grid[x_low][y_low][sliceNum+sliceStep] << "grid[x_low][y_high][sliceNum]=" << grid[x_low][y_high][sliceNum+sliceStep];
// //     //     qDebug() << "grid[x_high][y_low][sliceNum]=" << grid[x_high][y_low][sliceNum+sliceStep] << "grid[x_high][y_high][sliceNum]=" << grid[x_high][y_high][sliceNum+sliceStep];
// //     // }

// //     // Check bounds, if out of range, return default value
// //     if (x_low < 0 || x_high >= rows || y_low < 0 || y_high >= cols || x_low == x_high || y_low == y_high) {
// //         return default_value;
// //         qDebug() << "default-value is used!";
// //     }

// //     // if ((x_high >= rows) || (y_high >= cols)){
// //     //     x_ratio = 0;
// //     //     y_ratio = 0;

// //     // } else {

// //     double epsilon = 1e-10;
// //         // Bilinear interpolation
// //     if ( ((X[x_high][0] - X[x_low][0]) == 0) || (std::abs(X[x_high][0] - X[x_low][0]) < epsilon) ){
// //         x_ratio = 0;
// //     } else {
// //         x_ratio = (x - X[x_low][0]/mag) / (X[x_high][0]/mag - X[x_low][0]/mag);
// //     }
// //     //
// //     if ( ((Y[0][y_high] - Y[0][y_low]) == 0) || (std::abs(Y[0][y_high] - Y[0][y_low]) < epsilon) ){
// //         y_ratio = 0;
// //     } else {
// //         y_ratio = (y - Y[0][y_low]/mag) / (Y[0][y_high]/mag - Y[0][y_low]/mag);
// //     }
// //     // }

// //     double top = ( (grid[x_low][y_low][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_low][y_low][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_low][y_low][sliceNum+sliceStep]) ) * (1 - y_ratio)
// //                  + ( (grid[x_low][y_high][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_low][y_high][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_low][y_high][sliceNum+sliceStep]) ) * y_ratio;

// //     double bottom = ( (grid[x_high][y_low][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_high][y_low][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_high][y_low][sliceNum+sliceStep]) ) * (1 - y_ratio)
// //                     + ( (grid[x_high][y_high][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_high][y_high][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_high][y_high][sliceNum+sliceStep]) ) * y_ratio;

// //     return ((1 - x_ratio) * top + x_ratio * bottom); // still doesn't handle point 0,0, even with numInterp=1 : TODO

// // }

// // void ConeBeamSimulation::SubPixelInterpolation(const std::vector<std::vector<std::vector<int>>>& subject, std::vector<std::vector<double>>& beta,
// //                                     std::vector<std::vector<double>>& delta, const Materials::refractiveIndex& n1, const Materials::refractiveIndex& n2,
// //                                     const double& sourceToSubjectDist, const std::vector<double>& Mag, const int& sliceNum, const int& sliceStep, const int& m_numInterp)
// // {
// //     // double M_interp = 0;
// //     // beta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
// //     // delta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));

// //     // loop over different interpolation steps and accumulating them
// //     // for (int interp = 0 ; interp < m_numInterp ; interp++)
// //     // {
// //     // M_interp = m_SDD/(sourceToSubjectDist - thickness/2 + (sliceNum-1)*m_pixelSize/globalMag + interp*m_pixelSize/globalMag/m_numInterp);
// //     // M_interp = m_SDD/(sourceToSubjectDist + sliceNum*m_pixelSize/globalMag + interp*m_pixelSize/globalMag/m_numInterp);

// //     // M_interp = m_SDD/(sourceToSubjectDist + sliceNum*m_pixelSize/Mag[0] + (interp/m_numInterp)*(sliceStep)*m_pixelSize/Mag[0]);
// //     // qDebug() << "M_interp=" <<  M_interp;

// //     int numPixels = beta.size();

// //     for (int i = 0 ; i < numPixels ; i++){
// //         for (int j = 0 ; j < numPixels ; j++){
// //             beta[i][j] = ( (subject[i][j][sliceNum] == 1) ? ((n1.imaginaryPart) * subject[i][j][sliceNum]) : ((n2.imaginaryPart/2) * subject[i][j][sliceNum]) );
// //             delta[i][j] = ( (subject[i][j][sliceNum] == 1) ? ((n1.realPart) * subject[i][j][sliceNum]) : ((n2.realPart/2) * subject[i][j][sliceNum]) );
// //             // beta[i][j]+=Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.imaginaryPart, n2.imaginaryPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
// //             // delta[i][j]+=Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.realPart, n2.realPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
// //         }
// //     }
// //     // }

// //     // averaging over interpolations
// //     // for (int i = 0 ; i < m_numPixels ; i++){
// //     //     for (int j = 0 ; j < m_numPixels ; j++){
// //     //         beta[i][j] = beta[i][j] / m_numInterp;
// //     //         delta[i][j]= delta[i][j] / m_numInterp;
// //     //     }
// //     // }

// // }


// // void ConeBeamSimulation::PropagateInMaterial(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const std::vector<std::vector<std::vector<int>>>& subject,
// //                                     const double& sliceThickness, const double& sourceToSubjectDist, const std::vector<double>& Mag, const std::vector<double>& fresnelMag, const int& numMVSlices,
// //                                     const int& numVoxelSlicesInZ, const Materials::refractiveIndex& refrIndx1, const Materials::refractiveIndex& refrIndx2, const size_t& energyIndx)
// // {
// //     // declarations for each slice
// //     double coneBeamSliceThickness = 0, coneBeamSliceThickness2 = 0; // cone-beam slice thickness
// //     int sliceStep = numVoxelSlicesInZ/numMVSlices;

// //     // tempoarary loop matrices
// //     std::vector<std::vector<double>> beta; // imaginary part of the complex refractive index
// //     std::vector<std::vector<double>> delta; // real part of the complex refractive index
// //     std::vector<std::vector<double>> Gx_I; // intensity gradient along x : here vertical side
// //     std::vector<std::vector<double>> Gy_I; // intensity gradient along y : here horizontal side
// //     std::vector<std::vector<double>> Gx_phi; // phase gradient along x : here vertical side
// //     std::vector<std::vector<double>> Gy_phi; // phase gradient along y : here horizontal side
// //     std::vector<std::vector<double>> L_phi; // 2D Laplacian of phi

// //     // create magnified thickness mesh and interpolate refrative index in sub-pixel resolution
// //     for (int s = 0 ; s < numVoxelSlicesInZ ; s+=sliceStep)
// //     {
// //         beta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
// //         delta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
// //         Gx_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
// //         Gy_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
// //         Gx_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
// //         Gy_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
// //         L_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));

// //         SubPixelInterpolation(subject, beta, delta, refrIndx1, refrIndx2, sourceToSubjectDist, Mag, s, sliceStep, simPtr->m_numInterp); // sub-pixel interpolation of the complex refrative index
// //         Gradient2D(I, Gx_I, Gy_I, m_pixelSize/Mag[s+sliceStep]);
// //         Gradient2D(phi, Gx_phi, Gy_phi, m_pixelSize/Mag[s+sliceStep]);
// //         Laplacian2D(phi, L_phi, m_pixelSize/Mag[s+sliceStep]);

// //         // qDebug() << "m_waveNumber[energyIndx]" << m_waveNumber[energyIndx];
// //         // qDebug() << "sliceNum=" << s+sliceStep;
// //         // if (Mag[0] == 1){
// //         //     coneBeamSliceThickness = sliceThickness;
// //         //     for (int i = 0 ; i < m_numPixels ; i++){
// //         //         for (int j = 0 ; j < m_numPixels ; j++){
// //         //             I[i][j] = ( I[i][j] - (m_pixelSize/(simPtr->m_opticalMag)/(simPtr->m_waveNumber[energyIndx])) *
// //         //                                      (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) ) *
// //         //                       (std::exp(-2 * simPtr->m_waveNumber[energyIndx] * beta[i][j] * coneBeamSliceThickness));

// //         //             phi[i][j] = phi[i][j] - (simPtr->m_waveNumber[energyIndx] * delta[i][j] * coneBeamSliceThickness);
// //         //         }
// //         //     }
// //         // } else {
// //             // update cone-beam magnified slice thickness and then I and phi
// //             for (int i = 0 ; i < m_numPixels ; i++){
// //                 for (int j = 0 ; j < m_numPixels ; j++){
// //                     coneBeamSliceThickness = std::sqrt( simPtr->m_rsqr[i][j]/(std::pow(Mag[s+sliceStep], 2)) + std::pow((sourceToSubjectDist + (s+sliceStep)*m_pixelSize/Mag[0]), 2) )
// //                                              - std::sqrt( simPtr->m_rsqr[i][j]/(std::pow(Mag[s], 2)) + std::pow((sourceToSubjectDist + (s)*m_pixelSize/Mag[0]), 2) );
// //                     // OR
// //                     // coneBeamSliceThickness2 = std::sqrt(std::pow( (std::sqrt(m_rsqr[i][j])/(std::pow(Mag[s+sliceStep], 2))) - (std::sqrt(m_rsqr[i][j])/(std::pow(Mag[s], 2))) , 2) + (std::pow(sliceThickness*sliceStep, 2)));
// // //
// //                     I[i][j] = (1/fresnelMag[s+sliceStep])*( I[i][j] - (m_pixelSize/Mag[0]/fresnelMag[s+sliceStep]/(simPtr->m_waveNumber[energyIndx])) *
// //                                                                                (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) ) *
// //                               (std::exp(-2 * simPtr->m_waveNumber[energyIndx] * beta[i][j] * coneBeamSliceThickness));

// //                     phi[i][j] = phi[i][j] - (simPtr->m_waveNumber[energyIndx] * delta[i][j] * coneBeamSliceThickness);
// //                     if (coneBeamSliceThickness <= 0){
// //                         qDebug() << "******* CBSliceThickness is <= 0 at (" << i << "," << j << ")" << " *******";
// //                     }

// //                     // if ( I[i][j] < 0 ){
// //                     //     I[i][j] = std::abs((I[i-1][j] + I[i+1][j] + I[i][j-1] + I[i][j+1])/4);
// //                     //     // qDebug() << " ++++++ s=" << s << " ++++++";
// //                     //     // cc++;
// //                     //     // qDebug() << "I(" << i << "," << j << ")=" << I[i][j] << ", phi(" << i << "," << j << ")=" << phi[i][j];
// //                     //     // qDebug() << "coneBeamSliceThickness=" << coneBeamSliceThickness;
// //                     //     // qDebug() << "(i, j)=(" << i << "," << j << ") -> beta[i][j]=" << beta[i][j] << " ; delta[i][j]=" << delta[i][j];
// //                     //     // qDebug() << "GxI = " << Gx_I[i][j] << ", GyI = " << Gy_I[i][j] << ", Gxphi = " << Gx_phi[i][j] << ", Gyphi = " << Gy_phi[i][j] << ", Lphi = " << L_phi[j][j];
// //                     //     // qDebug() << "(I-a*b)*c, a->" << (m_pixelSize/Mag[0]/fresnelMag[s+sliceStep]/m_waveNumber[energyIndx]);
// //                     //     // qDebug() << "(I-a*b)*c, b->" << (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]);
// //                     //     // qDebug() << "(I-a*b)*c, c->" << std::exp(-2 * m_waveNumber[energyIndx] * beta[i][j] * coneBeamSliceThickness);
// //                     //     // qDebug() << "phi-d, d->" << m_waveNumber[energyIndx] * delta[i][j] * coneBeamSliceThickness;
// //                     //     // exit(1);
// //                     // }
// //                 }
// //             }
// //         // }
// //         qDebug() << "coneBeamSliceThickness=" << coneBeamSliceThickness << "; coneBeamSliceThickness2=" << coneBeamSliceThickness2;
// //     }

// // }

// // void ConeBeamSimulation::PropagateInFreeSpace(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const double& thickness1, const double& dist1, const double& dist2,
// //                                       const double& mag2, const double& mag1, const size_t& energyIndx, const bool isLastPropagation)
// // {
// //     // 1: first located element
// //     // thickness1 or 2: thickness of the first or second component from the source : i.e., object, diffuser, or detector(in this case=0)
// //     // dist1 or 2 : source to 1st or 2nd component distance : i.e. SOD, SDD, or SDD
// //     // Mag1 or 2 : relevant magnification factor(global/dynamic) of the 1st or 2nd located component (from the source)

// //     double energyWeight, detectorResponse = 0;
// //     double propagDist = 0;

// //     // tempoarary loop matrices
// //     std::vector<std::vector<double>> Gx_I(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // intensity gradient along x : here vertical side
// //     std::vector<std::vector<double>> Gy_I(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));  // intensity gradient along y : here horizontal side
// //     std::vector<std::vector<double>> Gx_phi(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // phase gradient along x : here vertical side
// //     std::vector<std::vector<double>> Gy_phi(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // phase gradient along y : here horizontal side
// //     std::vector<std::vector<double>> L_phi(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // 2D Laplacian of phi

// //     Gradient2D(I, Gx_I, Gy_I, m_pixelSize/mag1);
// //     Gradient2D(phi, Gx_phi, Gy_phi, m_pixelSize/mag1);
// //     Laplacian2D(phi, L_phi, m_pixelSize/mag1);

// //     if (isLastPropagation){
// //         energyWeight = simPtr->m_spectrumVector[energyIndx];
// //         detectorResponse = simPtr->m_detEnergyResponse[energyIndx];
// //         qDebug() << "It is the last propagation.\n\n";
// //     } else {
// //         energyWeight = 1;
// //         detectorResponse = 1;
// //     }

// //     // if ((mag1 == 1) && (mag2 == 1)){ // parallel beam

// //     //     propagDist = dist2 - (dist1 + thickness1);
// //     //     for (int i = 0 ; i < m_numPixels ; i++){
// //     //         for (int j = 0 ; j< m_numPixels ; j++){
// //     //             I[i][j] = energyWeight * detectorResponse * ( I[i][j] - (propagDist/m_waveNumber[energyIndx]) *
// //     //                                                                        (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) );
// //     //         }
// //     //     }

// //     // } else { // cone-beam
// //         for (int i = 0 ; i < m_numPixels ; i++){
// //             for (int j = 0 ; j< m_numPixels ; j++){

// //                 propagDist = std::sqrt(std::pow(dist2, 2) + (simPtr->m_rsqr[i][j])/(mag2*mag2))
// //                              - std::sqrt(std::pow((dist1+thickness1),2) + (simPtr->m_rsqr[i][j])/(mag1*mag1));

// //                 I[i][j] = energyWeight * detectorResponse * ( I[i][j] - (propagDist/dist2*dist1/(simPtr->m_waveNumber[energyIndx])) * (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) );

// //             }
// //         }
// //     // }
// // }

// void ConeBeamSimulation::Gradient2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& Gx, std::vector<std::vector<double>>& Gy, const double spacing)
// {
//     int rows = matrix.size();
//     int cols = matrix[0].size();

//     // Compute the gradient in the x-direction (here between rows)
//     for (int i = 1; i < rows - 1; ++i) {
//         for (int j = 0; j < cols; ++j) {
//             Gx[i][j] = (matrix[i + 1][j] - matrix[i - 1][j]) / (2.0 * spacing);
//         }
//     }
//     // Handle boundary conditions for Gx (forward and backward difference)
//     for (int j = 0; j < cols; ++j) {
//         Gx[0][j] = (matrix[1][j] - matrix[0][j]) / spacing;  // Forward difference for top row
//         Gx[rows - 1][j] = (matrix[rows - 1][j] - matrix[rows - 2][j]) / spacing;  // Backward difference for bottom row
//     }

//     // Compute the gradient in the y-direction (here between columns)
//     for (int i = 0; i < rows; ++i) {
//         for (int j = 1; j < cols - 1; ++j) {
//             Gy[i][j] = (matrix[i][j + 1] - matrix[i][j - 1]) / (2.0 * spacing);
//         }
//     }
//     // Handle boundary conditions for Gy (forward and backward difference)
//     for (int i = 0; i < rows; ++i) {
//         Gy[i][0] = (matrix[i][1] - matrix[i][0]) / spacing;  // Forward difference for left column
//         Gy[i][cols - 1] = (matrix[i][cols - 1] - matrix[i][cols - 2]) / spacing;  // Backward difference for right column
//     }

//     // avoid nan values: not at this stage
//     for (int i = 0 ; i < rows - 1 ; i++){
//         for (int j = 0 ; j < cols -1 ; j++){
//             if (std::isnan(Gx[i][j]) || std::isinf(Gx[i][j])){
//                 Gx[i][j] = 0;
//             }
//             if (std::isnan(Gy[i][j]) || std::isinf(Gy[i][j])){
//                 Gy[i][j] = 0;
//             }
//         }
//     }

// }

// void ConeBeamSimulation::Laplacian2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& L, const double& spacing)
// {
//     int rows = matrix.size();
//     int cols = matrix[0].size();

//     // Compute the Laplacian for the interior points (central difference)
//     for (int i = 1; i < rows - 1; ++i) {
//         for (int j = 1; j < cols - 1; ++j) {
//             // Second derivative in the x-direction
//             double d2x = (matrix[i + 1][j] - 2 * matrix[i][j] + matrix[i - 1][j]) / (spacing * spacing);

//             // Second derivative in the y-direction
//             double d2y = (matrix[i][j + 1] - 2 * matrix[i][j] + matrix[i][j - 1]) / (spacing * spacing);

//             // Laplacian is the sum of second derivatives
//             L[i][j] = d2x + d2y;
//         }
//     }

//     // Handle boundary conditions using forward/backward differences
//     // For the edges, we can use one-sided finite differences
//     // Top and bottom rows
//     for (int j = 0; j < cols; ++j) {
//         L[0][j] = (matrix[1][j] - 2 * matrix[0][j] + matrix[0][j]) / (spacing * spacing); // forward difference
//         L[rows - 1][j] = (matrix[rows - 1][j] - 2 * matrix[rows - 1][j] + matrix[rows - 2][j]) / (spacing * spacing); // backward difference
//     }

//     // Left and right columns
//     for (int i = 0; i < rows; ++i) {
//         L[i][0] = (matrix[i][1] - 2 * matrix[i][0] + matrix[i][0]) / (spacing * spacing); // forward difference
//         L[i][cols - 1] = (matrix[i][cols - 1] - 2 * matrix[i][cols - 1] + matrix[i][cols - 2]) / (spacing * spacing); // backward difference
//     }

//     // avoid nan values
//     for (int i = 0 ; i < rows - 1 ; i++){
//         for (int j = 0 ; j < cols -1 ; j++){
//             if (std::isnan(L[i][j]) || std::isinf(L[i][j])){
//                 L[i][j] = 0;
//             }
//         }
//     }


// }

//  std::vector<std::vector<std::vector<int>>> ConeBeamSimulation::rotate3DObjectY(std::vector<std::vector<std::vector<int>>>& object, const int& numProj)
// {
//     double rotAng = 360/numProj * pi /180.0; // radians

//     // rotation matrix around Y-axis
//     double cosTheta = std::cos(rotAng);
//     double sinTheta = std::sin(rotAng);

//     // size of 3D object
//     int Nx = object.size();
//     int Ny = object[0].size();
//     int Nz = object[0][0].size();

//     std::vector<std::vector<std::vector<int>>> rotatedObject(Nx, std::vector<std::vector<int>>(Ny, std::vector<int>(Nz, 0)));

//     // Get the center of the object (assume rotation around the center)
//     int cx = Nx / 2;
//     int cy = Ny / 2;
//     int cz = Nz / 2;

//     // Iterate through each point in the 3D rotated object -> backward mapping
//     for (int xr = 0; xr < Nx; ++xr) {
//         for (int yr = 0; yr < Ny; ++yr) {
//             for (int zr = 0; zr < Nz; ++zr) {

//                 // Translate point to the origin (center of object)
//                 int xt = xr - cx;
//                 int yt = yr - cy;
//                 int zt = zr - cz;

//                 // Apply the rotation matrix (only Y-axis rotation here)
//                 double xo = cosTheta * xt - sinTheta * zt;
//                 double yo = yt;  // No change in y for Y-axis rotation
//                 double zo = sinTheta * xt + cosTheta * zt;

//                 // Translate back to the original coordinate system
//                 int xo_int = std::round(xo + cx);
//                 int yo_int = std::round(yo + cy);
//                 int zo_int = std::round(zo + cz);

//                 // Check if the new coordinates are within bounds, if so, assign the value
//                 if (xo_int >= 0 && xo_int < Nx && yo_int >= 0 && yo_int < Ny && zo_int >= 0 && zo_int < Nz) {
//                     rotatedObject[xr][yr][zr] = object[xo_int][yo_int][zo_int];
//                 }
//             }
//         }
//     }
//     return rotatedObject;
// }

// double ConeBeamSimulation::Interpolation2D(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y,
//                                    const std::vector<std::vector<std::vector<int>>>& grid, const int& sliceNum, const int& sliceStep, const double& n1_part,
//                                    const double& n2_part, const double& mag, const double x, const double y, double default_value)
// {
//     int rows = grid.size();
//     int cols = grid[0].size();

//     double x_ratio = 0, y_ratio = 0;

//     // Find the four surrounding points in X and Y
//     int x_low = 0, x_high = 0, y_low = 0, y_high = 0;

//     // Find the correct x indices (x_low and x_high)
//     for (int i = 0; i < rows - 1; i++) {
//         if ( ((X[i][0]/mag) <= x) && ((X[i + 1][0]/mag) >= x) ) {
//             x_low = i;
//             x_high = i + 1;
//             break;
//         }
//     }
//     // Find the correct y indices (y_low and y_high)
//     for (int i = 0; i < cols - 1; i++) {
//         if ( ((Y[0][i]/mag) <= x) && ((Y[0][i + 1]/mag) >= x) ) {
//             y_low = i;
//             y_high = i + 1;
//             break;
//         }
//     }

//     // if ((grid[x_low][y_low][sliceNum+sliceStep] == 0) || (grid[x_low][y_low][sliceNum+sliceStep] == 0) || (grid[x_low][y_low][sliceNum+sliceStep] == 0) || (grid[x_low][y_low][sliceNum+sliceStep] == 0)){
//     //     qDebug() << "x=" << x << ";y=" << y << ";x_low=" << x_low << ";x_high=" << x_high << ";y_low=" << y_low << ";y_high=" << y_high << "";
//     //     qDebug() << "grid[x_low][y_low][sliceNum]=" << grid[x_low][y_low][sliceNum+sliceStep] << "grid[x_low][y_high][sliceNum]=" << grid[x_low][y_high][sliceNum+sliceStep];
//     //     qDebug() << "grid[x_high][y_low][sliceNum]=" << grid[x_high][y_low][sliceNum+sliceStep] << "grid[x_high][y_high][sliceNum]=" << grid[x_high][y_high][sliceNum+sliceStep];
//     // }

//     // Check bounds, if out of range, return default value
//     if (x_low < 0 || x_high >= rows || y_low < 0 || y_high >= cols || x_low == x_high || y_low == y_high) {
//         return default_value;
//         qDebug() << "default-value is used!";
//     }

//     // if ((x_high >= rows) || (y_high >= cols)){
//     //     x_ratio = 0;
//     //     y_ratio = 0;

//     // } else {

//     double epsilon = 1e-10;
//         // Bilinear interpolation
//     if ( ((X[x_high][0] - X[x_low][0]) == 0) || (std::abs(X[x_high][0] - X[x_low][0]) < epsilon) ){
//             x_ratio = 0;
//         } else {
//             x_ratio = (x - X[x_low][0]/mag) / (X[x_high][0]/mag - X[x_low][0]/mag);
//         }
//         //
//         if ( ((Y[0][y_high] - Y[0][y_low]) == 0) || (std::abs(Y[0][y_high] - Y[0][y_low]) < epsilon) ){
//             y_ratio = 0;
//         } else {
//             y_ratio = (y - Y[0][y_low]/mag) / (Y[0][y_high]/mag - Y[0][y_low]/mag);
//         }
//     // }

//     double top = ( (grid[x_low][y_low][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_low][y_low][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_low][y_low][sliceNum+sliceStep]) ) * (1 - y_ratio)
//                  + ( (grid[x_low][y_high][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_low][y_high][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_low][y_high][sliceNum+sliceStep]) ) * y_ratio;

//     double bottom = ( (grid[x_high][y_low][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_high][y_low][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_high][y_low][sliceNum+sliceStep]) ) * (1 - y_ratio)
//                     + ( (grid[x_high][y_high][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_high][y_high][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_high][y_high][sliceNum+sliceStep]) ) * y_ratio;

//         return ((1 - x_ratio) * top + x_ratio * bottom); // still doesn't handle point 0,0, even with numInterp=1 : TODO

// }

// void ConeBeamSimulation::SubPixelInterpolation(const std::vector<std::vector<std::vector<int>>>& subject, std::vector<std::vector<double>>& beta,
//                                        std::vector<std::vector<double>>& delta, const Materials::refractiveIndex& n1, const Materials::refractiveIndex& n2,
//                                        const double& sourceToSubjectDist, const std::vector<double>& Mag, const int& sliceNum, const int& sliceStep, const int& m_numInterp)
// {
//     // double M_interp = 0;
//     // beta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
//     // delta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));

//     // loop over different interpolation steps and accumulating them
//     // for (int interp = 0 ; interp < m_numInterp ; interp++)
//     // {
//         // M_interp = m_SDD/(sourceToSubjectDist - thickness/2 + (sliceNum-1)*m_pixelSize/globalMag + interp*m_pixelSize/globalMag/m_numInterp);
//         // M_interp = m_SDD/(sourceToSubjectDist + sliceNum*m_pixelSize/globalMag + interp*m_pixelSize/globalMag/m_numInterp);

//         // M_interp = m_SDD/(sourceToSubjectDist + sliceNum*m_pixelSize/Mag[0] + (interp/m_numInterp)*(sliceStep)*m_pixelSize/Mag[0]);
//         // qDebug() << "M_interp=" <<  M_interp;

//         for (int i = 0 ; i < m_numPixels ; i++){
//             for (int j = 0 ; j < m_numPixels ; j++){
//                 beta[i][j] = ( (subject[i][j][sliceNum] == 1) ? ((n1.imaginaryPart) * subject[i][j][sliceNum]) : ((n2.imaginaryPart/2) * subject[i][j][sliceNum]) );
//                 delta[i][j] = ( (subject[i][j][sliceNum] == 1) ? ((n1.realPart) * subject[i][j][sliceNum]) : ((n2.realPart/2) * subject[i][j][sliceNum]) );
//                 // beta[i][j]+=Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.imaginaryPart, n2.imaginaryPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
//                 // delta[i][j]+=Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.realPart, n2.realPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
//             }
//         }
//     // }

//     // averaging over interpolations
//     // for (int i = 0 ; i < m_numPixels ; i++){
//     //     for (int j = 0 ; j < m_numPixels ; j++){
//     //         beta[i][j] = beta[i][j] / m_numInterp;
//     //         delta[i][j]= delta[i][j] / m_numInterp;
//     //     }
//     // }

// }

// void ConeBeamSimulation::PropagateInMaterial(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const std::vector<std::vector<std::vector<int>>>& subject, const double& sliceThickness,
//                                      const double& sourceToSubjectDist, const std::vector<double>& Mag, const std::vector<double>& fresnelMag,const int& numMVSlices,
//                                      const int& numVoxelSlicesInZ, const Materials::refractiveIndex& refrIndx1, const Materials::refractiveIndex& refrIndx2, const size_t& energyIndx)
// {
//     // declarations for each slice
//     double coneBeamSliceThickness = 0, coneBeamSliceThickness2 = 0; // cone-beam slice thickness
//     int sliceStep = numVoxelSlicesInZ/numMVSlices;

//     // tempoarary loop matrices
//     std::vector<std::vector<double>> beta; // imaginary part of the complex refractive index
//     std::vector<std::vector<double>> delta; // real part of the complex refractive index
//     std::vector<std::vector<double>> Gx_I; // intensity gradient along x : here vertical side
//     std::vector<std::vector<double>> Gy_I; // intensity gradient along y : here horizontal side
//     std::vector<std::vector<double>> Gx_phi; // phase gradient along x : here vertical side
//     std::vector<std::vector<double>> Gy_phi; // phase gradient along y : here horizontal side
//     std::vector<std::vector<double>> L_phi; // 2D Laplacian of phi

//     // create magnified thickness mesh and interpolate refrative index in sub-pixel resolution
//     for (int s = 0 ; s < numVoxelSlicesInZ ; s+=sliceStep)
//     {
//         beta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
//         delta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
//         Gx_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
//         Gy_I.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
//         Gx_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
//         Gy_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
//         L_phi.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));

//         SubPixelInterpolation(subject, beta, delta, refrIndx1, refrIndx2, sourceToSubjectDist, Mag, s, sliceStep, m_numInterp); // sub-pixel interpolation of the complex refrative index
//         Gradient2D(I, Gx_I, Gy_I, m_pixelSize/Mag[s+sliceStep]);
//         Gradient2D(phi, Gx_phi, Gy_phi, m_pixelSize/Mag[s+sliceStep]);
//         Laplacian2D(phi, L_phi, m_pixelSize/Mag[s+sliceStep]);

//         // qDebug() << "m_waveNumber[energyIndx]" << m_waveNumber[energyIndx];
//         // qDebug() << "sliceNum=" << s+sliceStep;
//         if (Mag[0] == 1){
//             coneBeamSliceThickness = sliceThickness;
//             for (int i = 0 ; i < m_numPixels ; i++){
//                 for (int j = 0 ; j < m_numPixels ; j++){
//                     I[i][j] = ( I[i][j] - (m_pixelSize/(simPtr->m_opticalMag)/(simPtr->m_waveNumber[energyIndx])) *
//                                 (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) ) *
//                               (std::exp(-2 * (simPtr->m_waveNumber[energyIndx]) * beta[i][j] * coneBeamSliceThickness));

//                     phi[i][j] = phi[i][j] - ((simPtr->m_waveNumber[energyIndx]) * delta[i][j] * coneBeamSliceThickness);
//                 }
//             }
//         } else {
//         // update cone-beam magnified slice thickness and then I and phi
//             for (int i = 0 ; i < m_numPixels ; i++){
//                 for (int j = 0 ; j < m_numPixels ; j++){
//                     coneBeamSliceThickness = std::sqrt( (simPtr->m_rsqr[i][j])/(std::pow(Mag[s+sliceStep], 2)) + std::pow((sourceToSubjectDist + (s+sliceStep)*m_pixelSize/Mag[0]), 2) )
//                                              - std::sqrt( (simPtr->m_rsqr[i][j])/(std::pow(Mag[s], 2)) + std::pow((sourceToSubjectDist + (s)*m_pixelSize/Mag[0]), 2) );
//                     // OR
//                     // coneBeamSliceThickness2 = std::sqrt(std::pow( (std::sqrt(m_rsqr[i][j])/(std::pow(Mag[s+sliceStep], 2))) - (std::sqrt(m_rsqr[i][j])/(std::pow(Mag[s], 2))) , 2) + (std::pow(sliceThickness*sliceStep, 2)));

//                     I[i][j] = (1/fresnelMag[s+sliceStep])*( I[i][j] - (m_pixelSize/Mag[0]/fresnelMag[s+sliceStep]/(simPtr->m_waveNumber[energyIndx])) *
//                                              (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) ) *
//                               (std::exp(-2 * (simPtr->m_waveNumber[energyIndx]) * beta[i][j] * coneBeamSliceThickness));

//                     phi[i][j] = phi[i][j] - ((simPtr->m_waveNumber[energyIndx]) * delta[i][j] * coneBeamSliceThickness);
//                     if (coneBeamSliceThickness <= 0){
//                         qDebug() << "******* CBSliceThickness is <= 0 at (" << i << "," << j << ")" << " *******";
//                     }

//                     // if ( I[i][j] < 0 ){
//                     //     I[i][j] = std::abs((I[i-1][j] + I[i+1][j] + I[i][j-1] + I[i][j+1])/4);
//                     //     // qDebug() << " ++++++ s=" << s << " ++++++";
//                     //     // cc++;
//                     //     // qDebug() << "I(" << i << "," << j << ")=" << I[i][j] << ", phi(" << i << "," << j << ")=" << phi[i][j];
//                     //     // qDebug() << "coneBeamSliceThickness=" << coneBeamSliceThickness;
//                     //     // qDebug() << "(i, j)=(" << i << "," << j << ") -> beta[i][j]=" << beta[i][j] << " ; delta[i][j]=" << delta[i][j];
//                     //     // qDebug() << "GxI = " << Gx_I[i][j] << ", GyI = " << Gy_I[i][j] << ", Gxphi = " << Gx_phi[i][j] << ", Gyphi = " << Gy_phi[i][j] << ", Lphi = " << L_phi[j][j];
//                     //     // qDebug() << "(I-a*b)*c, a->" << (m_pixelSize/Mag[0]/fresnelMag[s+sliceStep]/m_waveNumber[energyIndx]);
//                     //     // qDebug() << "(I-a*b)*c, b->" << (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]);
//                     //     // qDebug() << "(I-a*b)*c, c->" << std::exp(-2 * m_waveNumber[energyIndx] * beta[i][j] * coneBeamSliceThickness);
//                     //     // qDebug() << "phi-d, d->" << m_waveNumber[energyIndx] * delta[i][j] * coneBeamSliceThickness;
//                     //     // exit(1);
//                     // }
//                 }
//             }
//         }
//         qDebug() << "coneBeamSliceThickness=" << coneBeamSliceThickness << "; coneBeamSliceThickness2=" << coneBeamSliceThickness2;
//     }

// }


// void ConeBeamSimulation::PropagateInFreeSpace(std::vector<std::vector<double>>& I, std::vector<std::vector<double>>& phi, const double& thickness1, const double& dist1, const double& dist2,
//                                       const double& mag2, const double& mag1, const size_t& energyIndx, const bool isLastPropagation)
// {
//     // 1: first located element
//     // thickness1 or 2: thickness of the first or second component from the source : i.e., object, diffuser, or detector(in this case=0)
//     // dist1 or 2 : source to 1st or 2nd component distance : i.e. SOD, SDD, or SDD
//     // Mag1 or 2 : relevant magnification factor(global/dynamic) of the 1st or 2nd located component (from the source)

//     double energyWeight, detectorResponse = 0;
//     double propagDist = 0;

//     // tempoarary loop matrices
//     std::vector<std::vector<double>> Gx_I(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // intensity gradient along x : here vertical side
//     std::vector<std::vector<double>> Gy_I(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));  // intensity gradient along y : here horizontal side
//     std::vector<std::vector<double>> Gx_phi(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // phase gradient along x : here vertical side
//     std::vector<std::vector<double>> Gy_phi(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // phase gradient along y : here horizontal side
//     std::vector<std::vector<double>> L_phi(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000)); // 2D Laplacian of phi

//     Gradient2D(I, Gx_I, Gy_I, m_pixelSize/mag1);
//     Gradient2D(phi, Gx_phi, Gy_phi, m_pixelSize/mag1);
//     Laplacian2D(phi, L_phi, m_pixelSize/mag1);

//     if (isLastPropagation){
//         energyWeight = simPtr->m_spectrumVector[energyIndx];
//         detectorResponse = simPtr->m_detEnergyResponse[energyIndx];
//         qDebug() << "It is the last propagation.\n\n";
//     } else {
//         energyWeight = 1;
//         detectorResponse = 1;
//     }

//     if ((mag1 == 1) && (mag2 == 1)){ // parallel beam

//         propagDist = dist2 - (dist1 + thickness1);
//         for (int i = 0 ; i < m_numPixels ; i++){
//             for (int j = 0 ; j< m_numPixels ; j++){
//                 I[i][j] = energyWeight * detectorResponse * ( I[i][j] - (propagDist/(simPtr->m_waveNumber[energyIndx])) *
//                             (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) );
//             }
//         }

//     } else { // cone-beam
//         for (int i = 0 ; i < m_numPixels ; i++){
//             for (int j = 0 ; j< m_numPixels ; j++){

//                 propagDist = std::sqrt(std::pow(dist2, 2) + (simPtr->m_rsqr[i][j])/(mag2*mag2))
//                              - std::sqrt(std::pow((dist1+thickness1),2) + (simPtr->m_rsqr[i][j])/(mag1*mag1));

//                 I[i][j] = energyWeight * detectorResponse * ( I[i][j] - (propagDist/dist2*dist1/(simPtr->m_waveNumber[energyIndx])) * (Gx_I[i][j] * Gx_phi[i][j] + Gy_I[i][j] * Gy_phi[i][j] + 4 * I[i][j] * L_phi[i][j]) );

//             }
//         }
//     }
// }

// void ConeBeamSimulation::InitiateSimulation()
// {

//     qDebug() << "pixelSize : " << m_pixelSize;
//     qDebug() << "numPixels : " << m_numPixels;
//     qDebug() << "obj-thickness : " << m_objThickness;
//     qDebug() << "numVoxelsInZObj : " << m_numObjVoxelsZ;
//     qDebug() << "numMVSlices_obj : " << m_numMVSlicesObj;
//     qDebug() << "diff-thickness : " << m_diffThickness;
//     qDebug() << "m_numMVSlicesDiff : " << m_numMVSlicesDiff;
//     qDebug() << "m_numDiffVoxelsZ : " << m_numDiffVoxelsZ;
//     qDebug() << "numInterp : " << simPtr->m_numInterp;
//     qDebug() << "numProj : " << simPtr->m_numProj;
//     qDebug() << "M_obj : " << simPtr->m_M_obj;
//     qDebug() << "M_diff : " << simPtr->m_M_diff;
//     qDebug() << "fM_obj : " << simPtr->m_fM_obj;
//     qDebug() << "fM_diff : " << simPtr->m_fM_diff;

//     int bg= 0, fg = 0;

//     // Images
//     std::vector<std::vector<float>> I_bg(m_numPixels, std::vector<float>(m_numPixels, 0.0));
//     std::vector<std::vector<float>> phi_bg(m_numPixels, std::vector<float>(m_numPixels, 0.0));
//     std::vector<std::vector<std::vector<float>>> I_fg(m_numPixels, std::vector<std::vector<float>>(m_numPixels, std::vector<float>(m_numProj, 0.0)));
//     std::vector<std::vector<std::vector<float>>> phi_fg(m_numPixels, std::vector<std::vector<float>>(m_numPixels, std::vector<float>(m_numProj, 0.0)));

//     // temporary matrices
//     std::vector<std::vector<double>> I_bg_t, phi_bg_t, I_fg_t, phi_fg_t;

//     // create object
//     std::unique_ptr<Object> object = std::make_unique<Object>(m_numObjVoxelsZ, m_numPixels, (m_pixelSize/simPtr->m_M_obj[0]), simPtr->tab_diffuser_object);
//     std::vector<std::vector<std::vector<int>>> obj(m_numPixels, std::vector<std::vector<int>>(m_numPixels, std::vector<int>(m_numObjVoxelsZ, 0)));
//     for (uint i = 0 ; i < simPtr->m_physList.size() ; i++)
//     {
//         qDebug() << "physlist : d = " << simPtr->m_physList[i].density << ", Z/A = " << simPtr->m_physList[i].Z_A << ", name = " << simPtr->m_physList[i].name;
//     }
//     std::vector<Materials::refractiveIndex> n_obj = simPtr->mat->RefractiveIndex(simPtr->m_energyVector, simPtr->m_physList.at(3).formula, simPtr->m_physList.at(3).density);

//     if (simPtr->tab_diffuser_object->m_virtualObject){
//         obj.clear();
//         obj.shrink_to_fit();
//         std::vector<std::vector<std::vector<float>>> vObj(m_numPixels, std::vector<std::vector<float>>(m_numPixels, std::vector<float>(m_numObjVoxelsZ)));
//         object->CreateObject(vObj);
//     } else {
//         object->CreateObject(obj);
//         // obj = object->rotate3DObjectZ(obj); // 45 deg rotated
//     }

//     // create diffuser
//     std::unique_ptr<Diffuser> diffuser = std::make_unique<Diffuser>(m_numPixels, m_numDiffVoxelsZ, m_pixelSize/(simPtr->m_M_diff[0]), simPtr->tab_diffuser_object);
//     std::vector<std::vector<std::vector<int>>> diff(m_numPixels, std::vector<std::vector<int>>(m_numPixels, std::vector<int>(m_numDiffVoxelsZ, 2)));
//     std::vector<Materials::refractiveIndex> n_diff = simPtr->mat->RefractiveIndex(simPtr->m_energyVector, simPtr->m_physList.at(1).formula, simPtr->m_physList.at(1).density);
//     std::vector<Materials::refractiveIndex> n_base = simPtr->mat->RefractiveIndex(simPtr->m_energyVector, simPtr->m_physList.at(2).formula, simPtr->m_physList.at(2).density);
//     diffuser->CreateDiffuser(diff); // create diffuser
//     qDebug() << "numGrits=" << diffuser->m_numGrits;


//     // ********************************************
//     // ****  Simulate image over projections  *****
//     // ********************************************
//     for (int proj = 0 ; proj < m_numProj ; proj++)
//     {
//         qDebug() << "========== proj = " << proj+1 << " ==========\n";
//         // iterate over all X-ray energies
//         for (size_t e = 0 ; e < simPtr->m_energyVector.size() ; e++)
//         {
//             qDebug() << "beta_obj=" << n_obj[e].imaginaryPart << "; delta_obj=" << n_obj[e].realPart << "\n";
//             qDebug() << "beta_diff=" << n_diff[e].imaginaryPart << "; delta_diff=" << n_diff[e].realPart << "\n";
//             qDebug() << "beta_base=" << n_base[e].imaginaryPart << "; delta_base=" << n_base[e].realPart << "\n";

//             // reset for the next energy bin
//             I_fg_t.assign(m_numPixels, std::vector<double>(m_numPixels, 1.0));
//             phi_fg_t.assign(m_numPixels, std::vector<double>(m_numPixels, 0.0));

//             // reference (background) image (only with diffuser) : only for one projection is calculated
//             if (proj == 0){
//                 I_bg_t.assign(m_numPixels, std::vector<double>(m_numPixels, 1.0));
//                 phi_bg_t.assign(m_numPixels, std::vector<double>(m_numPixels, 0.0));

//                 qDebug() << "Background image : propagating in the diffuser ...\n";
//                 PropagateInMaterial(I_bg_t, phi_bg_t, diff, m_diffThickness, m_SdD, simPtr->m_M_diff, simPtr->m_fM_diff, m_numMVSlicesDiff, m_numDiffVoxelsZ, n_diff[e], n_base[e], e); // propagate through the diffuser
//                 qDebug() << "Background image : propagating from the diffuser to the detector ...\n";
//                 PropagateInFreeSpace(I_bg_t, phi_bg_t, m_diffThickness, m_SdD, m_SDD, 1, simPtr->m_M_diff[m_numDiffVoxelsZ], e, true); // propagate from diffuser to detector
//             }

//             // foreground image (diffuser + object)
//             if (m_SOD < m_SdD){ // object is located before diffuser

//                 qDebug() << "Foreground image : propagating in the object ...\n";
//                 PropagateInMaterial(I_fg_t, phi_fg_t, obj, m_objThickness, m_SOD, simPtr->m_M_obj, simPtr->m_fM_obj, m_numMVSlicesObj, m_numObjVoxelsZ, n_obj[e], n_obj[e], e); // propagate through the object
//                 qDebug() << "Foreground image : propagating from the object to the diffuser ...\n";
//                 PropagateInFreeSpace(I_fg_t, phi_fg_t, m_objThickness, m_SOD, m_SdD, simPtr->m_M_diff[0], simPtr->m_M_obj[m_numObjVoxelsZ], e, false); // propagate from object to diffuser (or from diffuser to object)
//                 qDebug() << "Foreground image : propagating in the diffuser ...\n";
//                 PropagateInMaterial(I_fg_t, phi_fg_t, diff, m_diffThickness, m_SdD, simPtr->m_M_diff, simPtr->m_fM_diff, m_numMVSlicesDiff, m_numDiffVoxelsZ, n_diff[e], n_base[e], e); // propagate through the diffuser
//                 qDebug() << "Foreground image : propagating from the diffuser to the detector ...\n";
//                 PropagateInFreeSpace(I_fg_t, phi_fg_t, m_diffThickness, m_SdD, m_SDD, 1, simPtr->m_M_diff[m_numDiffVoxelsZ], e, true); // propagate from diffuser to detector

//             } else if (m_SOD > m_SdD) { // object is located after diffuser

//                 qDebug() << "Foreground image : propagating in the diffuser ...\n";
//                 PropagateInMaterial(I_fg_t, phi_fg_t, diff, m_diffThickness, m_SdD, simPtr->m_M_diff, simPtr->m_fM_diff, m_numMVSlicesDiff, m_numDiffVoxelsZ, n_diff[e], n_base[e], e); // propagate through the diffuser
//                 qDebug() << "Foreground image : propagating from the the diffuser to the object ...\n";
//                 PropagateInFreeSpace(I_fg_t, phi_fg_t, m_diffThickness, m_SdD, m_SOD, simPtr->m_M_obj[0], simPtr->m_M_diff[m_numDiffVoxelsZ], e, false); // propagate from diffuser to object
//                 qDebug() << "Foreground image : propagating in the object ...\n";
//                 PropagateInMaterial(I_fg_t, phi_fg_t, obj, m_objThickness, m_SOD, simPtr->m_M_obj, simPtr->m_fM_obj, m_numMVSlicesObj, m_numObjVoxelsZ, n_obj[e], n_obj[e], e); // propagate through the object
//                 qDebug() << "Foreground image : propagating from the object to the detector ...\n";
//                 PropagateInFreeSpace(I_fg_t, phi_fg_t, m_objThickness, m_SOD, m_SDD, 1, simPtr->m_M_obj[m_numObjVoxelsZ], e, true); // propagate from object to detector

//             }

//             // if image values are in a very low scale, it would result in 0.
//             // TODO : normalize all pixels to min or max and then convert it to float
//             if (proj == 0){ // diff image
//                 for (int i = 0 ; i < m_numPixels ; i++){
//                     for (int j = 0 ; j < m_numPixels ; j++){

//                         I_bg[i][j] += static_cast<float>(I_bg_t[i][j]);
//                         phi_bg[i][j] += static_cast<float>(phi_bg_t[i][j]);
//                         if (I_bg_t[i][j] < 0){
//                             bg++;
//                             // I_bg[i][j] = 0;
//                         }
//                     }
//                 }
//             }
//             // diff + obj image
//             for (int i = 0 ; i < m_numPixels ; i++){
//                 for (int j = 0 ; j < m_numPixels ; j++){

//                     I_fg[i][j][proj] += static_cast<float>(I_fg_t[i][j]);
//                     phi_fg[i][j][proj] += static_cast<float>(phi_fg_t[i][j]);
//                     if (I_fg_t[i][j] < 0){
//                         fg++;
//                         // I_fg[i][j] = 0;
//                     }
//                 }
//             }
//             qDebug() << "bg=" << bg << ", fg=" << fg << "\n";

//             bg = 0;
//             fg = 0;

//         } // end of energy loop

//         // Apply detector PSF to final images
//         // BlurImage();

//         // rotate object
//         obj = rotate3DObjectY(obj, m_numProj);

//     } // end of projections loop


//     // save background (bg) and foreground (fg) images to files
//     qDebug() << "writing images to files ... :";

//     // there is no rotation/tomo for bg images
//     std::string output_image_I_bg = "I_bg";
//     std::string output_image_phi_bg = "phi_bg";
//     simPtr->WriteToImage2D(I_bg, output_image_I_bg, 0);
//     simPtr->WriteToImage2D(phi_bg, output_image_phi_bg, 0);

//     qDebug() << "size of I_fg is :" << I_fg.size() << " x " << I_fg[0].size() << " x " << I_fg[0][0].size();

//     std::string output_image_I_fg = "I_fg";
//     std::string output_image_phi_fg = "phi_fg";
//     simPtr->WriteToImage3D(I_fg, output_image_I_fg);
//     simPtr->WriteToImage3D(phi_fg, output_image_phi_fg);

//     std::string output_image_obj = "object.nii";
//     std::string output_image_diff = "diffuser.nii";

//     std::vector<std::vector<std::vector<float>>> objj(m_numPixels, std::vector<std::vector<float>>(m_numPixels, std::vector<float>(m_numObjVoxelsZ)));
//     std::vector<std::vector<std::vector<float>>> difff(m_numPixels, std::vector<std::vector<float>>(m_numPixels, std::vector<float>(m_numDiffVoxelsZ)));


//     for (int i = 0 ; i < m_numPixels ; i++){
//         for (int j = 0 ; j < m_numPixels ; j++){
//             for (int k = 0 ; k < m_numObjVoxelsZ ; k++){
//                 objj[i][j][k] = static_cast<float>(obj[i][j][k]);
//             }
//         }
//     }

//     for (int i = 0 ; i < m_numPixels ; i++){
//         for (int j = 0 ; j < m_numPixels ; j++){
//             for (int k = 0 ; k < m_numDiffVoxelsZ ; k++){
//                 difff[i][j][k] = static_cast<float>(diff[i][j][k]);
//             }
//         }
//     }

//     simPtr->WriteToImage3D(objj, output_image_obj);
//     simPtr->WriteToImage3D(difff, output_image_diff);


// }
