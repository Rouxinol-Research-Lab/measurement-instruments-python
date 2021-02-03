#include "heterodynethread.h"
#include <qmath.h>
#include <random>
#include <toml.hpp>

double lorentzian(double A, double freq, double res_freq, double sigma)
{
	return (sigma / 2) / ((freq - res_freq)*(freq - res_freq) + (sigma / 2)*(sigma / 2)) / 3.1415;
}

void wave(double A, double freq, double fase, double* t, double* wavepoints, int N)
{

	for (int i = 0; i < N; i++)
	{
		wavepoints[i] = A * cos(freq*t[i] + fase);
	}

}

double VtodBm(double V)
{
	auto c = 10 * log10(20);
	return 20 * log10(V) + c;
}

void HeterodyneThread::simulate()
{

	double lower_bound = -0.0001;
	double upper_bound = 0.0001;
	std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
	std::uniform_real_distribution<double> unif_dut(lower_bound * 10, upper_bound * 100);
	std::default_random_engine re;

	

	int N = int(heterodyneSettings.sampleRate*heterodyneSettings.timeRange);
	double* t = new double[N];
	double* sourceI = new double[N];
	double* sourceQ = new double[N];
	double* dut = new double[N];

	for (int i = 0; i < N; i++)
	{
		{
			QMutexLocker locker(&m_mutex);
			if (m_stop) break;
		}

		t[i] = 0.1*i;
	}
	double I, Q;


	for (double freq = heterodyneSettings.startFrequency; freq < heterodyneSettings.stopFrequency; freq += (heterodyneSettings.stopFrequency- heterodyneSettings.startFrequency)/heterodyneSettings.nSteps )
	{
		I = 0;
		Q = 0;

		for (int i = 0; i < N; i++)
		{
			{
				QMutexLocker locker(&m_mutex);
				if (m_stop) break;
			}


			wave(0.001 + unif(re), freq, 0, t, sourceI, N);
			wave(0.001 + unif(re), freq, 3.1415 / 2, t, sourceQ, N);
			wave(lorentzian(1, freq, 3, 0.1) + unif_dut(re), freq, 0, t, dut, N);

		}


		for (int i = 0; i < N; i++)
		{
			{
				QMutexLocker locker(&m_mutex);
				if (m_stop) break;
			}

			I += sourceI[i] * dut[i] / N;
			Q += sourceQ[i] * dut[i] / N;
		}
			
		double result = VtodBm(4 * sqrt(I*I + Q * Q));
		emit signalDataPoint(freq, result);
	}

	delete t;
	delete sourceI;
	delete sourceQ;
	delete dut;
}

int capture_waveform(ViSession viOsc, char *buffer, double *points, int channel, int nPoints)
{
	int error;


	error = viPrintf(viOsc, ":WAV:FORMAT WORD\n");
	error = viPrintf(viOsc, ":WAV:SOUR CHAN%d, \n", channel);

	error = viPrintf(viOsc, "WAVEFORM:YINCREMENT?\n");
	error = viScanf(viOsc, "%t", buffer);
	qDebug("yincrement -> %s", buffer);

	std::string s = buffer;
	std::istringstream os(s);
	double yincrement;
	os >> yincrement;

	error = viPrintf(viOsc, "WAVEFORM:YORIGIN?\n");
	error = viScanf(viOsc, "%t", buffer);
	qDebug("yorigin -> %s", buffer);

	std::string s2 = buffer;
	std::istringstream os2(s2);
	double yorigin;
	os2 >> yorigin;


	error = viPrintf(viOsc, ":WAV:DATA? 1, \n");
	error = viScanf(viOsc, "%t", buffer);

	int numberHead = (int)(buffer[1] - 0x30);
	int numberOfBytes = 0;

	for (int i = 2; i < numberHead + 2; i++)
	{
		numberOfBytes += (int)(buffer[i] - 0x30) * pow(10, numberHead - 1 - i + 2);
	}

	qDebug("*buffer -> %s", buffer);
	qDebug("*buffer -> %c", buffer[0]);
	qDebug("*buffer -> %c", buffer[1]);



	char numberOfBytesPerPoint = 2;
	int n_points = numberOfBytes / numberOfBytesPerPoint;



	for (int i = 0; i < n_points; i++)
	{
		unsigned int byte1 = buffer[2 + numberHead + numberOfBytesPerPoint * i];
		unsigned int byte2 = buffer[2 + numberHead + 1 + numberOfBytesPerPoint * i];

		if (((0xff << 24) & byte2) >> 24 == 255)
			byte2 += 0x00000100;

		points[i] = ((signed int)((byte1 << 8) | byte2))*yincrement + yorigin;
	}

	return 0;
}

void HeterodyneThread::execute()
{
	qDebug("Thread id inside run %d", (int)QThread::currentThreadId());

	int error;
	
	ViSession session, viPSG1, viPSG2, viAttenuator, viOsc;
	char *buffer = new char[5000];

	qDebug("PSG1 address: %s", heterodyneSettings.Source1Address.c_str());
	qDebug("PSG2 address: %s", heterodyneSettings.Source2Address.c_str());
	qDebug("Attenuator address: %s", heterodyneSettings.AttenuatorAddress.c_str());
	qDebug("Oscilloscope address: %s", heterodyneSettings.OscilloscopeAddress.c_str());

	error = viOpenDefaultRM(&session);
	if (error != VI_SUCCESS)
	{
		qDebug("Error in locating resources");
	}

	error = viOpen(session, heterodyneSettings.Source1Address.c_str(), VI_NO_LOCK, 10000, &viPSG1);
	if (error != VI_SUCCESS)
	{
		qDebug("Error in opening PSG1");
	}

	error = viOpen(session, heterodyneSettings.Source2Address.c_str(), VI_NO_LOCK, 10000, &viPSG2);
	if (error != VI_SUCCESS)
	{
		qDebug("Error in opening PSG2");
	}

	error = viOpen(session, heterodyneSettings.AttenuatorAddress.c_str(), VI_NO_LOCK, 10000, &viAttenuator);
	if (error != VI_SUCCESS)
	{
		qDebug("Error in opening Attenuator");
	}

	error = viOpen(session, heterodyneSettings.OscilloscopeAddress.c_str(), VI_NO_LOCK, 10000, &viOsc);
	if (error != VI_SUCCESS)
	{
		qDebug("Error in opening Oscilloscope");
	}

	error = viPrintf(viPSG1, "*IDN?\n");
	error = viScanf(viPSG1, "%t", buffer);
	qDebug("*IDN? -> %s", buffer);


	error = viPrintf(viPSG2, "*IDN?\n");
	error = viScanf(viPSG2, "%t", buffer);
	qDebug("*IDN? -> %s", buffer);


	error = viPrintf(viAttenuator, "*IDN?\n");
	error = viScanf(viAttenuator, "%t", buffer);
	qDebug("*IDN? -> %s", buffer);


	error = viPrintf(viOsc, "*IDN?\n");
	error = viScanf(viOsc, "%t", buffer);
	qDebug("*IDN? -> %s", buffer);

	error = viPrintf(viPSG1, ":OUTP 0\n");
	error = viPrintf(viPSG2, ":OUTP 0\n");

	error = viPrintf(viPSG1, ":OUTP:MOD 0\n");
	error = viPrintf(viPSG2, ":OUTP:MOD 0\n");

	error = viPrintf(viPSG1, ":UNIT:POW DBM\n");
	error = viPrintf(viPSG2, ":UNIT:POW DBM\n");

	error = viPrintf(viPSG1, ":SOUR:POW:LEV:IMM:AMPL %d\n", heterodyneSettings.source1Amp);
	error = viPrintf(viPSG2, ":SOUR:POW:LEV:IMM:AMPL %d\n", heterodyneSettings.source2Amp);

	error = viPrintf(viAttenuator, ":ATT:BANK1:Y %d\n", heterodyneSettings.attenuation);

	error = viPrintf(viOsc, ":STOP\n");
	error = viPrintf(viOsc, ":TIM:RANG %f\n", heterodyneSettings.timeRange);
	error = viPrintf(viOsc, ":ACQ:SRAT:ANAL %f\n", heterodyneSettings.sampleRate);
	error = viPrintf(viOsc, ":SINGLE\n");
	sleep(5);
	error = viPrintf(viOsc, ":WAV:POIN?\n");
	error = viScanf(viOsc, "%t", buffer);


	std::string str = buffer;
	int nPoints = std::stoi(str);

	delete[] buffer;
	buffer = new char[2 * nPoints + 1000];


	double *Y1 = new double[nPoints];
	double *Y2 = new double[nPoints];
	double *Y3 = new double[nPoints];





	std::ofstream datafile(heterodyneSettings.filename);

	datafile << "Frequency" << "," << "I" << "," << "Q" << "\n";

	const double dfreq = (heterodyneSettings.stopFrequency - heterodyneSettings.startFrequency) / (heterodyneSettings.nSteps - 1);

	const double ifFreq = heterodyneSettings.ifFrequency;


	double I;
	double Q;

	error = viPrintf(viPSG1, ":OUTP 1\n");
	error = viPrintf(viPSG2, ":OUTP 1\n");

	for (double freq = heterodyneSettings.startFrequency; freq <= heterodyneSettings.stopFrequency; freq += dfreq)
	{


		I = 0;
		Q = 0;


		for (int j = 0; j < heterodyneSettings.averages; j++)
		{

			{
				QMutexLocker locker(&m_mutex);
				if (m_stop) break;
			}

			error = viPrintf(viPSG1, ":FREQ %f\n", freq + ifFreq);
			error = viPrintf(viPSG2, ":FREQ %f\n", freq);

			usleep(100000);

			error = viPrintf(viOsc, ":SINGLE\n");

			usleep(100000);



			capture_waveform(viOsc, buffer, Y1, heterodyneSettings.ChannelI, nPoints);
			capture_waveform(viOsc, buffer, Y2, heterodyneSettings.ChannelSignal, nPoints);
			capture_waveform(viOsc, buffer, Y3, heterodyneSettings.ChannelQ, nPoints);



			double aI = 0;
			double aQ = 0;
			for (int i = 0; i < nPoints; i++)
			{
				aI += Y1[i] * Y2[i] / nPoints;
				aQ += Y3[i] * Y2[i] / nPoints;
			}
			I += aI / heterodyneSettings.averages;
			Q += aQ / heterodyneSettings.averages;
		}


		double result = VtodBm(4 * sqrt(I*I + Q * Q));



		datafile << freq << "," << I << ',' << Q << "\n";




	}

	delete[] Y1;
	delete[] Y2;
	delete[] Y3;
	delete[] buffer;

	error = viPrintf(viPSG1, ":OUTPUT 0\n");
	error = viPrintf(viPSG2, ":OUTPUT 0\n");

	datafile.close();

	viClose(viPSG1);
	viClose(viPSG2);
	viClose(viAttenuator);
	viClose(viOsc);
	viClose(session);
}


void HeterodyneThread::run()
{
	if (heterodyneSettings.isSimulation)
	{
		emit signalLog("Simulation");
		simulate();
	}
	else
	{
		emit signalLog("Hardware");
		execute();
	}
	
}

void HeterodyneThread::stop()
{
	qDebug("Thread::stop called from main thread: %d", currentThreadId());
	QMutexLocker locker(&m_mutex);
	m_stop = true;
}

HeterodyneThread::HeterodyneThread(std::string dir)
{

	qDebug("Thread::stop called from main thread: %d", currentThreadId());
	qDebug(dir.c_str());


	const auto data = toml::parse(dir.c_str());
	
	const auto& instruments = toml::find(data, "settings","instruments");

	const auto toSimulate = toml::find<bool>(data, "settings", "simulation");
	heterodyneSettings.isSimulation = toSimulate;

	const auto filename = toml::find<std::string>(data, "settings", "savefile_name");
	heterodyneSettings.filename = filename;

	const auto att_address = toml::find<std::string>(instruments, "Attenuator");
	const auto source1_address = toml::find<std::string>(instruments, "Source1");
	const auto source2_address = toml::find<std::string>(instruments, "Source2");
	const auto osc_address = toml::find<std::string>(instruments, "Oscilloscope");
	heterodyneSettings.Source1Address = source1_address;
	heterodyneSettings.Source2Address = source2_address;
	heterodyneSettings.AttenuatorAddress = att_address;
	heterodyneSettings.OscilloscopeAddress = osc_address;

	const auto Channel_I_ref = toml::find<int>(instruments, "Channel_I_ref");
	heterodyneSettings.ChannelI = Channel_I_ref;

	const auto Channel_Q_ref = toml::find<int>(instruments, "Channel_Q_ref");
	heterodyneSettings.ChannelQ = Channel_Q_ref;

	const auto Channel_DUT_signal = toml::find<int>(instruments,"Channel_DUT_signal");
	heterodyneSettings.ChannelSignal = Channel_DUT_signal;


	const auto& experiment = toml::find(data, "experiment");

	const auto frequencies = toml::find<std::vector<double>>(experiment, "frequencies");
	heterodyneSettings.startFrequency = frequencies[0];
	heterodyneSettings.stopFrequency = frequencies[1];

	const auto frequency_nsteps = toml::find<int>(experiment, "frequency_nsteps");
	heterodyneSettings.nSteps = frequency_nsteps;

	const auto attenuation = toml::find<int>(experiment, "attenuation");
	heterodyneSettings.attenuation = attenuation;

	const auto averages = toml::find<int>(experiment, "averages");
	heterodyneSettings.averages = averages;

	const auto timeRange = toml::find<double>(experiment, "time_range");
	heterodyneSettings.timeRange = timeRange;

	const auto sampleRate = toml::find<double>(experiment, "sample_rate");
	heterodyneSettings.sampleRate = sampleRate;

	const auto ifFreq = toml::find<double>(experiment, "if_frequency");
	heterodyneSettings.ifFrequency = ifFreq;

	const auto source1Amp = toml::find<int>(experiment, "source1_amp");
	heterodyneSettings.source1Amp = source1Amp;

	const auto source2Amp = toml::find<int>(experiment, "source2_amp");
	heterodyneSettings.source2Amp = source2Amp;


}