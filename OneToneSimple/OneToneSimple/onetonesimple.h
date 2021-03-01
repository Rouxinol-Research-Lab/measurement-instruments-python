#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_onetonesimple.h"

class OneToneSimple : public QMainWindow
{
    Q_OBJECT

public:
    OneToneSimple(QWidget *parent = Q_NULLPTR);

public slots:
	void receivedDataPoint(double, double);
	void receivedLog(const char*);
	void defineRange(double, double);

private slots:
	void on_StopButton_clicked();
	

private:
    Ui::OneToneSimpleClass ui;
	QVector<double> x, y;

signals:
	void stopMeasurement();
};
