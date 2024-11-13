#include "Diffuser.h"

#include <vector>
#include <random>
#include <algorithm>
#include <cmath>


Diffuser::Diffuser(const int numPixels, const int numDiffVoxelsZ, const double pixelDiff, DiffuserAndObject *diffuserAndObject)
    : m_numPixels(numPixels)
    , m_numDiffVoxelsZ(numDiffVoxelsZ)
    , m_pixelDiff(pixelDiff)
    , m_diffuserAndObject(diffuserAndObject)
{
    m_diffuserSize = {m_numPixels, m_numPixels, m_numDiffVoxelsZ};
    m_numGrits = CalculateNumGrits();

}

Diffuser::~Diffuser(){

}

int Diffuser::CalculateNumGrits()
{ // the volume of all average-size grit spheres to the base volume is 30-70%
    double baseVol = (m_numPixels*m_pixelDiff) * (m_numPixels*m_pixelDiff) * (m_numDiffVoxelsZ*m_pixelDiff);
    double gritVol = 4/3*pi*std::pow(((m_diffuserAndObject->getGritSize())/2), 3);

    if (m_diffuserAndObject->getGritDensity() =="Dense") {
        return (round(baseVol/gritVol*0.5));

    } else if (m_diffuserAndObject->getGritDensity() == "Standard") {
        return (round(baseVol/gritVol*0.2));

    } else if (m_diffuserAndObject->getGritDensity() == "Sparse") {
        return (round(baseVol/gritVol*0.1));
    } else {
        return 1;
    }
}

// generate random integer numbers for radius, center_x,y,z of the grit volume
void Diffuser::getRandomValue(std::vector<int>& myVec, int lower_size, int upper_size)
{
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> gritDis(lower_size, upper_size);

    for (int i = 0 ; i < m_numGrits ; i++)
    {
        myVec.push_back(gritDis(gen));
        // std::cout << "c: " << myVec[i] << "\n";
    }
}

void Diffuser::DistributeGrits(std::vector<int>& cx, std::vector<int>& cy, std::vector<int>& cz, std::vector<int>& r)
{
    int gritRadiusIndx = 0;
    cx.clear(); cy.clear(); cz.clear(); r.clear();

    double meanGritRadius = (m_diffuserAndObject->getGritSize()); // mean
    double stdDev = 0.2 * meanGritRadius; // std
    double maxGritRadius = meanGritRadius + 3 * stdDev;
    int maxGritRadiusIndx = maxGritRadius/m_pixelDiff;
    maxGritRadiusIndx = std::min(maxGritRadiusIndx, m_numDiffVoxelsZ);

    qDebug() << "Mean radius of the grits : " << meanGritRadius;
    qDebug() << "Radius stdDev of the grits  : " << stdDev;
    qDebug() << "Maximum grit radius : " << maxGritRadius;
    qDebug() << "Maximum grit radius index : " << maxGritRadiusIndx;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(meanGritRadius, stdDev);

    // Generate the radii
    for (int i = 0; i < m_numGrits; i++) {

        gritRadiusIndx = std::round((dis(gen))/m_pixelDiff);
        gritRadiusIndx = std::clamp(gritRadiusIndx, 1, maxGritRadiusIndx);
        r.push_back(gritRadiusIndx);

        // std::cout << "r= " << gritRadiusIndx << "\n";
    }

    getRandomValue(cx, maxGritRadiusIndx-1, std::max(m_diffuserSize[0]-maxGritRadiusIndx, maxGritRadiusIndx-1));
    getRandomValue(cy, maxGritRadiusIndx-1, std::max(m_diffuserSize[1]-maxGritRadiusIndx, maxGritRadiusIndx-1));
    getRandomValue(cz, maxGritRadiusIndx-1, std::max(m_diffuserSize[2]-maxGritRadiusIndx, maxGritRadiusIndx-1));
    // getRandomValue(cx, 0, m_diffuserSize[0]-1);
    // getRandomValue(cy, 0, m_diffuserSize[1]-1);
    // getRandomValue(cz, 0, m_diffuserSize[2]-1);
}

void Diffuser::CreateDiffuser(std::vector<std::vector<std::vector<int>>>& diffuser)
{
    std::vector<int> cx{0}, cy{0}, cz{0}, r{0};
    int r_squared = 0;

    DistributeGrits(cx, cy, cz, r);

    for (size_t s = 0; s < cx.size(); s++) {
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
                        diffuser[i][j][k] = 1; // for grit
                    // } else {
                    //     diffuser[i][j][k] = 2; // for base
                    }
                }
            }
        }
    }
}
