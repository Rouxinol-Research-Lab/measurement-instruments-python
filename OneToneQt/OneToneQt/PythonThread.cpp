#include "PythonThread.h"


double lorentzian(double A, double freq, double res_freq, double sigma)
{
	return (sigma / 2) / ((freq - res_freq)*(freq - res_freq) + (sigma / 2)*(sigma / 2))/3.1415;
}


void PythonThread::run()
{
	std::ofstream datafile(PythonSettings.filename);

	datafile << "---------------" << "Low Power" << "---------------" << "\n";

	datafile << "Frequency" << "," << "Amplitude" << "\n";

	double dfreq = (PythonSettings.stopFrequency - PythonSettings.startFrequency) / (PythonSettings.nSteps-1);
	qDebug("Thread id inside run %d", (int)QThread::currentThreadId());

	double half_freq = (PythonSettings.stopFrequency + PythonSettings.startFrequency) / 2;


	//low power
	for (double freq = PythonSettings.startFrequency; freq <= PythonSettings.stopFrequency; freq+= dfreq)
	{
		{
			QMutexLocker locker(&m_mutex);
			if (m_stop) break;
		}

		double result = lorentzian(10,freq, half_freq - 0.05f, 0.01f);

		usleep(10000);

		datafile << freq << "," << result << "\n";

		emit signalDataPoint(0,freq, result);
	}


	datafile << "---------------" << "High Power" << "---------------" << "\n";
	datafile << "Frequency" << "," << "Amplitude" << "\n";

	//high power
	for (double freq = PythonSettings.startFrequency; freq <= PythonSettings.stopFrequency; freq += dfreq)
	{
		{
			QMutexLocker locker(&m_mutex);
			if (m_stop) break;
		}

		double result = lorentzian(10, freq, half_freq, 0.01f);

		usleep(10000);

		datafile << freq << "," << result << "\n";

		emit signalDataPoint(1, freq, result);
	}

	datafile.close();

}

void PythonThread::loadAndStart(MeasurementSetting settings)
{
	{
		QMutexLocker locker(&m_mutex);
		m_stop = false;
	}
	PythonSettings = settings;
	emit this->start();
}


void PythonThread::stop()
{
	qDebug("Thread::stop called from main thread: %d", currentThreadId());
	QMutexLocker locker(&m_mutex);
	m_stop = true;
}