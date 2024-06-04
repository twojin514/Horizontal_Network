#include "Horizontal_Network.h"
#include "ui_Horizontal_Network.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QThread>

Horizontal_Network::Horizontal_Network(QWidget* parent)
    : QMainWindow(parent)
    , ui(new Ui::Horizontal_NetworkClass)
{
    ui->setupUi(this);

    // 시그널과 슬롯 연결
    connect(ui->controlButton, &QPushButton::clicked, this, &Horizontal_Network::on_controlButton_clicked);
    connect(ui->measureButton, &QPushButton::clicked, this, &Horizontal_Network::on_measureButton_clicked);
    connect(ui->angleButton, &QPushButton::clicked, this, &Horizontal_Network::on_angleButton_clicked);
    connect(ui->lineButton, &QPushButton::clicked, this, &Horizontal_Network::on_lineButton_clicked);
    connect(ui->startButton, &QPushButton::clicked, this, &Horizontal_Network::on_startButton_clicked);

}

Horizontal_Network::~Horizontal_Network()
{
    delete ui;
}

QString Horizontal_Network::getControlFilePath() const {
    return ui->ControlLineEdit->text();
}

QString Horizontal_Network::getMeasureFilePath() const {
    return ui->MeasureLineEdit->text();
}

QString Horizontal_Network::getAngleFilePath() const {
    return ui->AngleLineEdit->text();
}

QString Horizontal_Network::getLineFilePath() const {
    return ui->LineLineEdit->text();
}

QString Horizontal_Network::getTitle() const {
    return ui->titleLineEdit->text();
}

QString Horizontal_Network::getBusiness() const {
    return ui->businessLineEdit->text();
}

QString Horizontal_Network::getName() const {
    return ui->nameLineEdit->text();
}

QString Horizontal_Network::getOutput() const {
    return ui->outputLineEdit->text();
}


void Horizontal_Network::on_controlButton_clicked()
{
    QString fileName = getControlFilePath();

    if (!fileName.isEmpty()) {
        loadCSV(fileName, ui->ControlTable);
    }
}


void Horizontal_Network::on_measureButton_clicked()
{
    QString fileName = ui->MeasureLineEdit->text();

    if (!fileName.isEmpty()) {
        loadCSV(fileName, ui->MeasureTable);
    }
}

void Horizontal_Network::on_angleButton_clicked()
{
    QString fileName = ui->AngleLineEdit->text();

    if (!fileName.isEmpty()) {
        loadCSV(fileName, ui->AngleTable);
    }
}

void Horizontal_Network::on_lineButton_clicked()
{
    QString fileName = ui->LineLineEdit->text();

    if (!fileName.isEmpty()) {
        loadCSV(fileName, ui->LineTable);
    }
}

void Horizontal_Network::loadCSV(const QString& fileName, QTableWidget* tableWidget)
{
    QFile file(fileName);
    if (!file.exists()) {
        QMessageBox::warning(this, "File Error", "File does not exist: " + fileName);
        return;
    }

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::warning(this, "File Error", "Unable to open file: " + fileName);
        return;
    }

    QTextStream in(&file);
    QStringList headers = in.readLine().split(",");
    if (headers.isEmpty()) {
        QMessageBox::warning(this, "File Error", "Failed to read headers from file: " + fileName);
        return;
    }

    tableWidget->setColumnCount(headers.size());
    tableWidget->setHorizontalHeaderLabels(headers);

    tableWidget->setRowCount(0);
    int row = 0;
    while (!in.atEnd()) {
        QStringList fields = in.readLine().split(",");
        if (fields.size() != headers.size()) {
            QMessageBox::warning(this, "File Error", "Data format error in file: " + fileName);
            return;
        }
        tableWidget->insertRow(row);
        for (int col = 0; col < fields.size(); ++col) {
            tableWidget->setItem(row, col, new QTableWidgetItem(fields.at(col)));
        }
        ++row;
    }

    file.close();

}


void Horizontal_Network::on_startButton_clicked()
{
    // 프로그레스 바 초기화 및 상태 텍스트 업데이트
    ui->progressBar->setValue(0);
    ui->statusLabel->setText("Start.");

    emit startButtonClicked();

    // 프로그램이 진행되는 동안 프로그레스 바 업데이트
    for (int i = 0; i <= 100; ++i) {
        QCoreApplication::processEvents(); // 이벤트 루프 처리
        ui->progressBar->setValue(i);
        QThread::msleep(10); // 진행을 확인하기 위해 잠시 대기
    }

    // 작업이 완료되면 프로그레스 바 및 상태 텍스트 업데이트
    ui->progressBar->setValue(100);
    ui->statusLabel->setText("End.");
}

void Horizontal_Network::updateProgress(int value)
{
    ui->progressBar->setValue(value);
}

void Horizontal_Network::updateStatus(const QString& status)
{
    ui->statusLabel->setText(status);
}