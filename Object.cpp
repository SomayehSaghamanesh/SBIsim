#include "Object.h"

#include <QDebug>


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

void Object::CreateObject(std::vector<std::vector<std::vector<int>>>& object)
{
    if (m_diffuserAndObject->getObjectName() == "Cylinder"){
        CreateCylinder(object);

    } else if (m_diffuserAndObject->getObjectName() == "Sphere"){
        CreateSphere(object);

    } else {
        std::runtime_error("ERROR: Please select a valid object.");
        return;
    }
}

void Object::CreateObject(std::vector<std::vector<std::vector<float>>>& object)
{
    double mag_vObject = m_diffuserAndObject->getVObjectMag();
    object = m_diffuserAndObject->getVirtualObject(mag_vObject, m_numPixels, m_numObjVoxelsZ, m_pixelObj);
}

void Object::CreateSphere(std::vector<std::vector<std::vector<int>>>& object)
{
    int r = 0;

    if (m_diffuserAndObject){
        r = (m_diffuserAndObject->getSphereDiameter()/2)/m_pixelObj; // sphere radius as index
    }
    // double d = getSphereDiameter();
    // int r = (d/2)/m_pixelObj; // sphere radius as index
    int originInPlane = m_numPixels/2;
    int originZ = m_numObjVoxelsZ/2;

    qDebug() << "Sphere diameter: " << (m_diffuserAndObject->getSphereDiameter());
    // qDebug() << "m_pixepObj inside Object : " << m_pixelObj;
    qDebug() << "Sphere radius index : " << r;
    qDebug() << "Sphere origin in x-y : " << originInPlane;

    for (int i = 0 ; i < m_numPixels ; i++)
    {
        for (int j = 0 ; j < m_numPixels ; j++)
        {
            for (int  k = 0 ; k < m_numObjVoxelsZ ; k++)
            {
                if ( ((i - originInPlane)*(i - originInPlane) + (j - originInPlane)*(j - originInPlane) + (k - originZ)*(k - originZ)) <= r*r){

                    object[i][j][k] = 1;
                }

            }
        }
    }
}

void Object::CreateCylinder(std::vector<std::vector<std::vector<int>>>& object)
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

                    object[i][j][k] = 1;
                }

            }
        }
    }
}

// void Object::CreateVirtualObject(std::vector<std::vector<std::vector<float>>>& object)
// {
//     double mag_vObject = m_diffuserAndObject->getVObjectMag();
//     object = m_diffuserAndObject->getVirtualObject(mag_vObject, m_numPixels, m_numObjVoxelsZ, m_pixelObj);
// }
std::vector<std::vector<std::vector<int>>> Object::rotate3DObjectZ(std::vector<std::vector<std::vector<int>>>& object)
{
    const double pi = M_PI;
    double rotAng = 90 * pi /180.0; // radians

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

    // Iterate through each point in the 3D rotated object
    for (int xr = 0; xr < Nx; ++xr) {
        for (int yr = 0; yr < Ny; ++yr) {
            for (int zr = 0; zr < Nz; ++zr) {

                // Translate point back to the origin (center of object)
                int xt = xr - cx;
                int yt = yr - cy;
                int zt = zr - cz;

                // Apply the rotation matrix (only Z-axis rotation here)
                double xo = cosTheta * xt + sinTheta * yt;
                double yo = -sinTheta * xt + cosTheta * yt;
                double zo = zt; // No change in z for Z-axis rotation

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
