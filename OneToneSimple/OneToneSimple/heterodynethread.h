#pragma once

#include <QtCore>
#include <QMutex>
#include "onetonedatamodel.h"
#include <fstream>
#include "visa.h"


class HeterodyneThread : public QThread
{
	Q_OBJECT


private:
	QMutex m_mutex;
	bool m_stop = false;
	MeasurementSetting heterodyneSettings;
	void simulate();
	void execute();

protected:
	virtual void run();

public:
	HeterodyneThread(std::string);

public slots:
	void stop();

signals:
	void signalDataPoint(double, double);
	void signalLog(const char*);
};
