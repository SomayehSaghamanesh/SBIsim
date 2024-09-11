#include "Diffuser.h"(

#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

Diffuser::Diffuser()
{
    m_diffuserSize = {m_numPixels, m_numPixels, m_diffThicknessAsIndex};
    m_numGrit = 100;
}

Diffuser::~Diffuser(){}


// generate random integer numbers for radius, center_x,y,z of the grit volume
void Diffuser::getRandomValue(std::vector<int>& myVec, int lower_size, int upper_size)
{
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> gritDis(lower_size, upper_size);

    for (int i = 0 ; i < m_numGrit ; i++)
    {
        myVec.push_back(gritDis(gen));
    }
}

void Diffuser::DistributeGrits(std::vector<int>& cx, std::vector<int>& cy, std::vector<int>& cz, std::vector<int>& r)
{
    m_diffuserSize = {m_numPixels, m_numPixels, m_diffThicknessAsIndex};

    float meanGritRadius = getGritSizeMM()/2; // mean
    float stdDev = 0.2 * meanGritRadius; // std
    float maxGritRadius = meanGritRadius + 3 * stdDev;
    int maxGritRadiusIndx = maxGritRadius/m_pixelDiff;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(meanGritRadius, stdDev);
    // Generate the radii
    for (int i = 0; i < m_numGrit; i++) {
        int gritRadiusIndx = std::clamp(static_cast<int>(round(dis(gen))), 1, maxGritRadiusIndx);
        r.push_back(gritRadiusIndx);
    }

    getRandomValue(cx, maxGritRadiusIndx, m_diffuserSize[0]-maxGritRadiusIndx);
    getRandomValue(cy, maxGritRadiusIndx, m_diffuserSize[1]-maxGritRadiusIndx);
    getRandomValue(cz, maxGritRadiusIndx, m_diffuserSize[2]-maxGritRadiusIndx);
}

void Diffuser::CreateDiffuser(std::vector<std::vector<std::vector<double>>>& diffuser, double refrIndx)
{
    std::vector<int> cx, cy, cz, r;

    DistributeGrits(cx, cy, cz, r);

    for (size_t s = 0; s < cx.size(); ++s) {
        int r_squared = r[s] * r[s]; // Calculate squared radius

        // bounding box for the current sphere
        int cx_min = std::max(0, cx[s] - r[s]);
        int cx_max = std::min(m_diffuserSize[0] - 1, cx[s] + r[s]);
        int cy_min = std::max(0, cy[s] - r[s]);
        int cy_max = std::min(m_diffuserSize[1] - 1, cy[s] + r[s]);
        int cz_min = std::max(0, cz[s] - r[s]);
        int cz_max = std::min(m_diffuserSize[2] - 1, cz[s] + r[s]);

        // Iterate within the bounding box
        for (int i = cx_min; i <= cx_max; ++i) {
            for (int j = cy_min; j <= cy_max; ++j) {
                for (int k = cz_min; k <= cz_max; ++k) {
                    int dx = i - cx[s];
                    int dy = j - cy[s];
                    int dz = k - cz[s];
                    if (dx * dx + dy * dy + dz * dz < r_squared) {
                        diffuser[i][j][k] = refrIndx;
                    }
                }
            }
        }
    }
}

//     int Diffuser::getRandomValue(double weight1, double weight2, std::mt19937 &gen)
// {
//     std::uniform_real_distribution<> dis(0.0, 1.0);
//     return (dis(gen) < weight1 / (weight1 + weight2)) ? 1 : 2;
// }

// void Diffuser::initializeGritClusters(std::vector<std::vector<std::vector<int>>> &grid, int sizeX, int sizeY, int sizeZ, int avgDiameter)
// {
//     std::mt19937 gen(std::random_device{}());
//     std::uniform_int_distribution<> disX(0, sizeX - 1);
//     std::uniform_int_distribution<> disY(0, sizeY - 1);
//     std::uniform_int_distribution<> disZ(0, sizeZ - 1);

//     int numSeeds = (sizeX * sizeY * sizeZ) / (avgDiameter * avgDiameter * avgDiameter);
//     std::vector<std::tuple<int, int, int>> seedPositions;

//     for (int i = 0; i < numSeeds; ++i) {
//         int x = disX(gen);
//         int y = disY(gen);
//         int z = disZ(gen);

//         int value = getRandomValue(1.0, 1.0, gen); // Initially equal probability for 1 or 2
//         grid[x][y][z] = value;
//         seedPositions.push_back({x, y, z});
//     }

//     std::queue<std::tuple<int, int, int>> queue;
//     for (const auto &pos : seedPositions) {
//         queue.push(pos);
//     }

//     while (!queue.empty()) {
//         auto [x, y, z] = queue.front();
//         queue.pop();
//         int currentValue = grid[x][y][z];

//         std::vector<std::tuple<int, int, int>> neighbors = {
//             {x - 1, y, z}, {x + 1, y, z}, {x, y - 1, z}, {x, y + 1, z},
//             {x, y, z - 1}, {x, y, z + 1}
//         };

//         for (const auto &[nx, ny, nz] : neighbors) {
//             if (nx >= 0 && ny >= 0 && nz >= 0 && nx < sizeX && ny < sizeY && nz < sizeZ && grid[nx][ny][nz] == 0) {
//                 double weightForSame = 2.0; // Higher probability to be the same as the current cluster
//                 double weightForOther = 1.0;
//                 int newValue = getRandomValue(weightForSame, weightForOther, gen);

//                 grid[nx][ny][nz] = newValue == 1 ? currentValue : (3 - currentValue); // Ensure only 1 or 2 are used
//                 queue.push({nx, ny, nz});
//             }
//         }
//     }
// }

// void Diffuser::VisualizeDiffuser()
// {
//     int sizeX = 20, sizeY = 20, sizeZ = 20;
//     int avgDiameter = 5; // Target average diameter of clusters

//     // Step 1: Initialize a 3D vector
//     std::vector<std::vector<std::vector<int>>> grid(sizeX, std::vector<std::vector<int>>(sizeY, std::vector<int>(sizeZ, 0)));

//     // Step 2: Initialize and grow clusters
//     initializeGritClusters(grid, sizeX, sizeY, sizeZ, avgDiameter);

//     // Output the generated 3D vector (for demonstration purposes, we'll print a 2D slice)
//     for (int z = 0; z < sizeZ; ++z) {
//         std::cout << "Slice " << z << ":\n";
//         for (int y = 0; y < sizeY; ++y) {
//             for (int x = 0; x < sizeX; ++x) {
//                 std::cout << grid[x][y][z] << " ";
//             }
//             std::cout << std::endl;
//         }
//         std::cout << std::endl;
//     }

// }
