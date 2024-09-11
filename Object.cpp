#include "Object.h"


Object::Object() {

}

Object::~Object()
{

}


void Object::CreateSphere(std::vector<std::vector<std::vector<double>>>& object, double refrIndx)
{
    int numPixelsZ = getObjectThickness()/m_pixelObj;
    int r = round(getSphereDiameter()/2);

    int originInPlane = floor(m_numPixels/2);

    for (int i = 0 ; i < m_numPixels ; i++)
    {
        for (int j = 0 ; j < m_numPixels ; j++)
        {
            for (int  k = 0 ; k < numPixelsZ ; k++)
            {
                if ( ((i - originInPlane)*(i - originInPlane) + (j - originInPlane)*(j - originInPlane) + (k - originInPlane)*(k - originInPlane)) <= r*r){

                    object[i][j][k] = refrIndx;
                }

            }
        }
    }
}

void Object::CreateCylinder(std::vector<std::vector<std::vector<double>>>& object, double refrIndx)
{
    float d = getCylinderDiameter();
    int r = d/2;
    int h = getCylinderHeight()/2; // half-height
    int numPixelsZ = d/m_pixelObj;
    int originZ = floor(numPixelsZ/2);
    int originInPlane = floor(m_numPixels/2);

    for (int i = 0 ; i < m_numPixels ; i++)
    {
        for (int j = 0 ; j < m_numPixels ; j++)
        {
            for (int  k = 0 ; k < numPixelsZ ; k++)
            {
                if ( (((i - originInPlane)*(i - originInPlane) + (k - originZ)*(k - originZ)) <= r*r) && (std::abs(j - originInPlane) <= h/2) ){

                    object[i][j][k] = refrIndx;
                }

            }
        }
    }
}

