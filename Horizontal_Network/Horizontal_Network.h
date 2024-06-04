#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtWidgets/QMainWindow>
#include "ui_Horizontal_Network.h"
#include <QString>
#include <QThread>
#include <QTableWidget>

QT_BEGIN_NAMESPACE
namespace Ui { class Horizontal_NetworkClass; };
QT_END_NAMESPACE

class Horizontal_Network : public QMainWindow
{
    Q_OBJECT

public:
    Horizontal_Network(QWidget *parent = nullptr);
    ~Horizontal_Network();
    QString getControlFilePath() const;
    QString getMeasureFilePath() const;
    QString getAngleFilePath() const;
    QString getLineFilePath() const;
    QString getTitle() const;
    QString getBusiness() const;
    QString getName() const;
    QString getOutput() const;


signals:
    void startButtonClicked();

private slots:
    void on_controlButton_clicked();
    void on_measureButton_clicked();
    void on_angleButton_clicked();
    void on_lineButton_clicked();
    void on_startButton_clicked();
    void updateProgress(int value);
    void updateStatus(const QString& status);

private:
    Ui::Horizontal_NetworkClass* ui;
    void loadCSV(const QString& fileName, QTableWidget* tableWidget);
    QString controlFilePath;
    QString measureFilePath;
    QString angleFilePath;
    QString lineFilePath;
    QString title;
    QString business;
    QString name;
    QString output;
    QThread workerThread;
};

#endif // HORIZONTAL_NETWORKCLASS_H