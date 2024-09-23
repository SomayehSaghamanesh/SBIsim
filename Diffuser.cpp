#include "Diffuser.h"
#include "Materials.h"

#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

Diffuser::Diffuser(int numPixels, int numDiffVoxelsZ, double pixelDiff, DiffuserAndObject *diffuserAndObject)
    : m_numPixels(numPixels)
    , m_numDiffVoxelsZ(numDiffVoxelsZ)
    , m_pixelDiff(pixelDiff)
    , m_diffuserAndObject(diffuserAndObject)
{
    m_diffuserSize = {m_numPixels, m_numPixels, m_numDiffVoxelsZ};
    m_numGrit = 100;
}

Diffuser::~Diffuser(){

}


// generate random integer numbers for radius, center_x,y,z of the grit volume
void Diffuser::getRandomValue(std::vector<int>& myVec, int lower_size, int upper_size)
{
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> gritDis(lower_size, upper_size);

    for (int i = 0 ; i < m_numGrit ; i++)
    {
        myVec.push_back(gritDis(gen));
        qDebug() << "c: " << myVec[i];
    }
}

void Diffuser::DistributeGrits(std::vector<int>& cx, std::vector<int>& cy, std::vector<int>& cz, std::vector<int>& r)
{
    int gritRadiusIndx = 0;
    cx.clear(); cy.clear(); cz.clear(); r.clear();

    double meanGritRadius = (m_diffuserAndObject->getGritSizeMM()); // mean
    double stdDev = 0.2 * meanGritRadius; // std
    double maxGritRadius = meanGritRadius + 3 * stdDev;
    int maxGritRadiusIndx = maxGritRadius/m_pixelDiff;
    maxGritRadiusIndx = std::min(maxGritRadiusIndx, m_numDiffVoxelsZ);

    qDebug() << "meanGritRadius : " << meanGritRadius;
    qDebug() << "stdDev : " << stdDev;
    qDebug() << "maxGritRadius : " << maxGritRadius;
    qDebug() << "maxGritRadiusIndx : " << maxGritRadiusIndx;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(meanGritRadius, stdDev);

    // Generate the radii
    for (int i = 0; i < m_numGrit; i++) {

        gritRadiusIndx = std::round((dis(gen))/m_pixelDiff);
        gritRadiusIndx = std::clamp(gritRadiusIndx, 1, maxGritRadiusIndx);
        r.push_back(gritRadiusIndx);

        // qDebug() << "r= " << gritRadiusIndx;
    }

    getRandomValue(cx, maxGritRadiusIndx, m_diffuserSize[0]-maxGritRadiusIndx);
    getRandomValue(cy, maxGritRadiusIndx, m_diffuserSize[1]-maxGritRadiusIndx);
    getRandomValue(cz, maxGritRadiusIndx, m_diffuserSize[2]-maxGritRadiusIndx);
}

void Diffuser::CreateDiffuser(std::vector<std::vector<std::vector<Materials::refractiveIndex>>>& diffuser, Materials::refractiveIndex& n_diff, Materials::refractiveIndex& n_base)
{
    std::vector<int> cx{0}, cy{0}, cz{0}, r{0};
    int r_squared = 0;

    DistributeGrits(cx, cy, cz, r);

    for (size_t s = 0; s < cx.size(); ++s) {
        r_squared = r[s] * r[s]; // Calculate squared radius

        // bounding box for the current sphere
        int cx_min = std::max(0, cx[s] - r[s]); // to avoid negative index
        int cx_max = std::min(m_diffuserSize[0] - 1, cx[s] + r[s]); // to avoid outlier index
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
                    if (dx * dx + dy * dy + dz * dz <= r_squared) {
                        diffuser[i][j][k] = n_diff;
                    } else {
                        diffuser[i][j][k] = n_base;
                    }
                }
            }
        }
    }
}
