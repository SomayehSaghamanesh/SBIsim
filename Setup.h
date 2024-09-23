#ifndef SETUP_H
#define SETUP_H

#include <QWidget>

namespace Ui {
class Setup;
}

class Setup : public QWidget
{
    Q_OBJECT

public:
    explicit Setup(QWidget *parent = nullptr);
    ~Setup();

    // double getSourceToObjectDist();
    // double getSourcetoDiffuserDist();
    // double getSourcetoDetectorDist();
    void getDistances();
    int getNumProjections();
    double getTomoAngleDeg();

    bool m_radioFlag{false}, m_tomoFlag{false};
    double m_SOD{0}, m_SDD{0}, m_SdD{0};

public slots :


private slots:

    void on_radioButton_Radiography_toggled(bool checked);
    void on_radioButton_Tomography_toggled(bool checked);


private:
signals:
    void StartButtonClicked();


private:
    Ui::Setup *ui;

};

#endif // SETUP_H
