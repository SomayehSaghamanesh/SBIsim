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
    double getPixelSizeMM();
    double getDetectorThickness();
    QString getDetectorMaterial();
    int getNumPixels();
    double getFOVmm();
    std::vector<double> getWaveNumber(std::vector<double>& energy);
    void Meshgrid(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y, int& numPixels, double& pixelSize);
    void DetectorCoordinates(std::vector<std::vector<double>>& X, std::vector<std::vector<double>>& Y, std::vector<std::vector<double>>& rsqr, int& numPixels);


    bool m_parBeam{false};
    bool m_coneBeam{false};
    bool m_monochrome{false};
    bool m_polychrome{false};

    // double m_pixelSizeMM{0};
    // int m_numPixels{0};

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
