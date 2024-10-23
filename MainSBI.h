#ifndef MAINSBI_H
#define MAINSBI_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui {
class MainSBI;
}
QT_END_NAMESPACE

class MainSBI : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainSBI(QWidget *parent = nullptr);
    ~MainSBI();

private slots:
    void on_actionSimulation_triggered();

    void on_actionExit_triggered();

private:
    Ui::MainSBI *ui;
};
#endif // MAINSBI_H
