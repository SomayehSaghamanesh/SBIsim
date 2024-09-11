#ifndef DIFFUSERANDOBJECT_H
#define DIFFUSERANDOBJECT_H

#include <QWidget>
#include <QVector>

namespace Ui {
class DiffuserAndObject;
}

class DiffuserAndObject : public QWidget
{
    Q_OBJECT

    public:
        explicit DiffuserAndObject(QWidget *parent = nullptr);
        ~DiffuserAndObject();

        float getCylinderDiameter();
        float getCylinderHeight();
        float getSphereDiameter();
        QVector<QVector<QVector<float>>> getVirtualObject();
        QString getObjectMaterial();
        int getNumMVSlices();
        int getNumInterpolations();
        float getGritSizeMM();
        QString getGritMaterial();
        QString getBaseMaterial();
        double getObjectThickness();
        double getDiffuserThickness();
        void MagFactors();


    public slots:


    private slots:

        void on_radioButton_Cylinder_toggled(bool checked);
        void on_radioButton_Sphere_toggled(bool checked);
        void on_radioButton_VObject_toggled(bool checked);
        void on_pushButton_VObject_clicked();


    private:

        Ui::DiffuserAndObject *ui;

        bool m_cylinder{false}, m_sphere{false}, m_virtualObject{false};
        QVector<QVector<QVector<float>>> m_vObject;
        QStringList m_materialsList = {};
        QStringList m_num_MVS_slices = {"few", "medium", "large"};


    protected :

        double m_M_obj{1}, m_M_diff{1};
        std::vector<double> m_dM_obj{1}, m_dM_diff{1}, m_fM_obj{1}, m_fM_diff{1};
        double m_objThickness{0}, m_diffThickness{0}, m_pixelSize{0}, m_pixelObj{0}, m_pixelDiff{0};
        int m_objThicknessAsIndex{0}, m_diffThicknessAsIndex{0}, m_numPixels{0};

        // float m_gritSize;
        // QString m_objectMaterial{""}, m_gritMaterial{""}, m_baseMaterial{""};

        void MaterialsList();

    };


#endif // DIFFUSERANDOBJECT_H