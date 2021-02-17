#include "onetonesimple.h"
#include <QtWidgets/QApplication>
#include "heterodynethread.h"


int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	OneToneSimple w;
	w.show();

	bool ok;
	const char* filename = "C:\\Users\\lucam\\Documents\\Mestrado\\Medidas\\2021-02-17 Analuse 2 do circuito heterodino\\test.toml";
	QString dirText = QInputDialog::getText(0, "Experiment TOML File",
		"File name:", QLineEdit::Normal,
		filename, &ok);
	
	if(ok) {

		qDebug("Thread id %d", (int)QThread::currentThreadId());
		
		HeterodyneThread t(dirText.toUtf8().toStdString());
		QObject::connect(&t, SIGNAL(signalDataPoint(double, double)), &w, SLOT(receivedDataPoint(double, double)));
		QObject::connect(&t, SIGNAL(signalLog(char*)), &w, SLOT(receivedLog(char*)));
		QObject::connect(&w, SIGNAL(stopMeasurement()), &t, SLOT(stop()));

		t.start();
		return a.exec();
	}

	return 0;
}
