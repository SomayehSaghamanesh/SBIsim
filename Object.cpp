#include "Object.h"


Object::Object(int numObjVoxelsZ, int numPixels, double pixelObj, DiffuserAndObject* diffuserAndObject)
    : m_numObjVoxelsZ(numObjVoxelsZ)
    , m_numPixels(numPixels)
    , m_pixelObj(pixelObj)
    , m_diffuserAndObject(diffuserAndObject)
{

}

Object::~Object()
{
}

void Object::CreateSphere(std::vector<std::vector<std::vector<Materials::refractiveIndex>>>& object, Materials::refractiveIndex& n_obj)
{
    int r = 0;

    if (m_diffuserAndObject){
        r = (m_diffuserAndObject->getSphereDiameter()/2)/m_pixelObj; // sphere radius as index
    }
    // double d = getSphereDiameter();
    // int r = (d/2)/m_pixelObj; // sphere radius as index
    int originInPlane = m_numPixels/2;

    qDebug() << "sphere diam : " << (m_diffuserAndObject->getSphereDiameter());
    qDebug() << "m_pixepObj inside Object : " << m_pixelObj;
    qDebug() << "sphere radius index : " << r;
    qDebug() << "sphere originInPlane : " << originInPlane;

    for (int i = 0 ; i < m_numPixels ; i++)
    {
        for (int j = 0 ; j < m_numPixels ; j++)
        {
            for (int  k = 0 ; k < m_numObjVoxelsZ ; k++)
            {
                if ( ((i - originInPlane)*(i - originInPlane) + (j - originInPlane)*(j - originInPlane) + (k - originInPlane)*(k - originInPlane)) <= r*r){

                    object[i][j][k] = n_obj;
                }

            }
        }
    }
}

void Object::CreateCylinder(std::vector<std::vector<std::vector<Materials::refractiveIndex>>>& object, Materials::refractiveIndex& n_obj)
{
    int r = (m_diffuserAndObject->getCylinderDiameter())/2/m_pixelObj;
    int h = (m_diffuserAndObject->getCylinderHeight())/2/m_pixelObj; // half-height
    int originZ = m_numObjVoxelsZ/2;
    int originInPlane = m_numPixels/2;

    for (int i = 0 ; i < m_numPixels ; i++)
    {
        for (int j = 0 ; j < m_numPixels ; j++)
        {
            for (int  k = 0 ; k < m_numObjVoxelsZ ; k++)
            {
                if ( (((i - originInPlane)*(i - originInPlane) + (k - originZ)*(k - originZ)) <= r*r) && (std::abs(j - originInPlane) <= h) ){

                    object[i][j][k] = n_obj;
                }

            }
        }
    }
}

