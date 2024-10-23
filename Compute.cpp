#include "Compute.h"

Compute::Compute() {}

Compute::~Compute() {}


void Compute::Gradient2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& Gx, std::vector<std::vector<double>>& Gy, const double spacing)
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

void Compute::Laplacian2D(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& L, const double& spacing)
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

std::vector<std::vector<std::vector<int>>> Compute::rotate3DObjectY(std::vector<std::vector<std::vector<int>>>& object, const int& numProjs)
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

    // Iterate through each point in the 3D rotated object -> backward mapping
    for (int xr = 0; xr < Nx; ++xr) {
        for (int yr = 0; yr < Ny; ++yr) {
            for (int zr = 0; zr < Nz; ++zr) {

                // Translate point to the origin (center of object)
                int xt = xr - cx;
                int yt = yr - cy;
                int zt = zr - cz;

                // Apply the rotation matrix (only Y-axis rotation here)
                double xo = cosTheta * xt - sinTheta * zt;
                double yo = yt;  // No change in y for Y-axis rotation
                double zo = sinTheta * xt + cosTheta * zt;

                // Translate back to the original coordinate system
                int xo_int = std::round(xo + cx);
                int yo_int = std::round(yo + cy);
                int zo_int = std::round(zo + cz);

                // Check if the new coordinates are within bounds, if so, assign the value
                if (xo_int >= 0 && xo_int < Nx && yo_int >= 0 && yo_int < Ny && zo_int >= 0 && zo_int < Nz) {
                    rotatedObject[xr][yr][zr] = object[xo_int][yo_int][zo_int];
                }
            }
        }
    }
    return rotatedObject;
}

double Compute::Interpolation2D(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y,
                                   const std::vector<std::vector<std::vector<int>>>& grid, const int& sliceNum, const int& sliceStep, const double& n1_part,
                                   const double& n2_part, const double& mag, const double x, const double y, double default_value)
{
    int rows = grid.size();
    int cols = grid[0].size();

    double x_ratio = 0, y_ratio = 0;

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

    // if ((x_high >= rows) || (y_high >= cols)){
    //     x_ratio = 0;
    //     y_ratio = 0;

    // } else {

    double epsilon = 1e-10;
        // Bilinear interpolation
    if ( ((X[x_high][0] - X[x_low][0]) == 0) || (std::abs(X[x_high][0] - X[x_low][0]) < epsilon) ){
        x_ratio = 0;
    } else {
        x_ratio = (x - X[x_low][0]/mag) / (X[x_high][0]/mag - X[x_low][0]/mag);
    }
    //
    if ( ((Y[0][y_high] - Y[0][y_low]) == 0) || (std::abs(Y[0][y_high] - Y[0][y_low]) < epsilon) ){
        y_ratio = 0;
    } else {
        y_ratio = (y - Y[0][y_low]/mag) / (Y[0][y_high]/mag - Y[0][y_low]/mag);
    }
    // }

    double top = ( (grid[x_low][y_low][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_low][y_low][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_low][y_low][sliceNum+sliceStep]) ) * (1 - y_ratio)
                 + ( (grid[x_low][y_high][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_low][y_high][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_low][y_high][sliceNum+sliceStep]) ) * y_ratio;

    double bottom = ( (grid[x_high][y_low][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_high][y_low][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_high][y_low][sliceNum+sliceStep]) ) * (1 - y_ratio)
                    + ( (grid[x_high][y_high][sliceNum+sliceStep] == 1) ? ((n1_part) * grid[x_high][y_high][sliceNum+sliceStep]) : ((n2_part/2) * grid[x_high][y_high][sliceNum+sliceStep]) ) * y_ratio;

    return ((1 - x_ratio) * top + x_ratio * bottom); // still doesn't handle point 0,0, even with numInterp=1 : TODO

}

void Compute::SubPixelInterpolation(const std::vector<std::vector<std::vector<int>>>& subject, std::vector<std::vector<double>>& beta,
                                       std::vector<std::vector<double>>& delta, const Materials::refractiveIndex& n1, const Materials::refractiveIndex& n2,
                                       const double& sourceToSubjectDist, const std::vector<double>& Mag, const int& sliceNum, const int& sliceStep, const int& m_numInterp)
{
    // double M_interp = 0;
    // beta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));
    // delta.assign(m_numPixels, std::vector<double>(m_numPixels, 0.000000000000000));

    // loop over different interpolation steps and accumulating them
    // for (int interp = 0 ; interp < m_numInterp ; interp++)
    // {
    // M_interp = m_SDD/(sourceToSubjectDist - thickness/2 + (sliceNum-1)*m_pixelSize/globalMag + interp*m_pixelSize/globalMag/m_numInterp);
    // M_interp = m_SDD/(sourceToSubjectDist + sliceNum*m_pixelSize/globalMag + interp*m_pixelSize/globalMag/m_numInterp);

    // M_interp = m_SDD/(sourceToSubjectDist + sliceNum*m_pixelSize/Mag[0] + (interp/m_numInterp)*(sliceStep)*m_pixelSize/Mag[0]);
    // qDebug() << "M_interp=" <<  M_interp;

    int numPixels = beta.size();

    for (int i = 0 ; i < numPixels ; i++){
        for (int j = 0 ; j < numPixels ; j++){
            beta[i][j] = ( (subject[i][j][sliceNum] == 1) ? ((n1.imaginaryPart) * subject[i][j][sliceNum]) : ((n2.imaginaryPart/2) * subject[i][j][sliceNum]) );
            delta[i][j] = ( (subject[i][j][sliceNum] == 1) ? ((n1.realPart) * subject[i][j][sliceNum]) : ((n2.realPart/2) * subject[i][j][sliceNum]) );
            // beta[i][j]+=Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.imaginaryPart, n2.imaginaryPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
            // delta[i][j]+=Interpolation2D(m_X, m_Y, subject, sliceNum, sliceStep, n1.realPart, n2.realPart, Mag[sliceNum], m_X[i][j]/M_interp, m_Y[i][j]/M_interp, 0);
        }
    }
    // }

    // averaging over interpolations
    // for (int i = 0 ; i < m_numPixels ; i++){
    //     for (int j = 0 ; j < m_numPixels ; j++){
    //         beta[i][j] = beta[i][j] / m_numInterp;
    //         delta[i][j]= delta[i][j] / m_numInterp;
    //     }
    // }

}
