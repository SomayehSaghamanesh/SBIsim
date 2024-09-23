#ifndef SOURCEANDDETECTOR_H
#define SOURCEANDDETECTOR_H

#include <QWidget>

namespace Ui {
class SourceAndDetector;
}

class SourceAndDetector : public QWidget
{
    Q_OBJECT

public:
    explicit SourceAndDetector(QWidget *parent = nullptr);
    ~SourceAndDetector();

    void getXrayEnergy(std::vector<double>& energyVector, std::vector<double>& spectrumVector);
    void getPixelSizeMM();
    double getDetectorThickness();
    QString getDetectorMaterial();
    int getNumPixels();
    double getFOVmm();
    double waveLength(double energy);

    bool m_parBeam{false};
    bool m_coneBeam{false};
    bool m_monochrome{false};
    bool m_polychrome{false};

    double m_pixelSizeMM{0};

    std::vector<double> m_waveNumber();

public slots:


private slots:

    void on_pushButton_energySpectrum_clicked();

    void on_radioButton_parBeam_toggled(bool checked);

    void on_radioButton_coneBeam_toggled(bool checked);

    void on_radioButton_Monochrom_toggled(bool checked);

    void on_radioButton_Polychrom_toggled(bool checked);



private:
    Ui::SourceAndDetector *ui;

    // double m_monoEnergy = 0.0;
    // QVector<double> m_energyVector{0};
    // QVector<double> m_SpectrumVector{0};
    QStringList m_materialsList = {};
    // QString m_detectorMaterial = "";
    // double m_detectorThickness = 0.0;
    // double m_pixelSize = 0.0;

    const double hc = 1.24e-9; // [keV.m]
    const double pi = M_PI;


    void MaterialsList();


};

#endif // SOURCEANDDETECTOR_H
