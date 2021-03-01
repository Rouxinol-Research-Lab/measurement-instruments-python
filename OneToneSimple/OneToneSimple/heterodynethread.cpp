#include "heterodynethread.h"
#include <qmath.h>
#include <toml.hpp>
#include <string>

double lorentzian(double A, double freq, double res_freq, double sigma)
{
	return (sigma / 2) / ((freq - res_freq)*(freq - res_freq) + (sigma / 2)*(sigma / 2)) / 3.1415;
}

double VtodBm(double V)
{
	auto c = 10 * log10(20);
	return 20 * log10(V) + c;
}

void HeterodyneThread::simulate()
{
	emit signalLog("Simulation");


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


			sourceI[i] =  cos(freq*t[i]);
			sourceQ[i] = cos(freq*t[i]+ 3.1415 / 2 );
			dut[i] = lorentzian(1, freq, 5.6e9, 1e6)*cos(freq*t[i]);

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
	int error = 0;


	error = viPrintf(viOsc, ":WAV:FORMAT WORD\n");
	if (error != VI_SUCCESS)
	{
		return error;
	}
	error = viPrintf(viOsc, ":WAV:SOUR CHAN%d, \n", channel);
	if (error != VI_SUCCESS)
	{
		return error;
	}

	error = viPrintf(viOsc, "WAVEFORM:YINCREMENT?\n");
	if (error != VI_SUCCESS)
	{
		return error;
	}
	error = viScanf(viOsc, "%t", buffer);
	if (error != VI_SUCCESS)
	{
		return error;
	}

	std::string s = buffer;
	std::istringstream os(s);
	double yincrement;
	os >> yincrement;

	error = viPrintf(viOsc, "WAVEFORM:YORIGIN?\n");
	if (error != VI_SUCCESS)
	{
		return error;
	}
	error = viScanf(viOsc, "%t", buffer);
	if (error != VI_SUCCESS)
	{
		return error;
	}

	std::string s2 = buffer;
	std::istringstream os2(s2);
	double yorigin;
	os2 >> yorigin;


	error = viPrintf(viOsc, ":WAV:DATA? 1, \n");
	if (error != VI_SUCCESS)
	{
		return error;
	}
	error = viScanf(viOsc, "%t", buffer);
	if (error != VI_SUCCESS)
	{
		return error;
	}

	int numberHead = (int)(buffer[1] - 0x30);
	int numberOfBytes = 0;

	for (int i = 2; i < numberHead + 2; i++)
	{
		numberOfBytes += (int)(buffer[i] - 0x30) * pow(10, numberHead - 1 - i + 2);
	}
/*
	qDebug("*buffer -> %s", buffer);
	qDebug("*buffer -> %c", buffer[0]);
	qDebug("*buffer -> %c", buffer[1]);
	*/


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

	return error;
}

void HeterodyneThread::execute()
{
	emit signalLog("Instruments.");

	int error;
	
	ViSession session, viPSG1, viPSG2, viAttenuator, viOsc;
	char *buffer = new char[5000];

	std::string Address_Source1 = "PSG1 address: ";
	Address_Source1 += heterodyneSettings.Source1Address.c_str();
	emit signalLog(Address_Source1.c_str());

	std::string Address_Source2 = "PSG2 address: ";
	Address_Source2 += heterodyneSettings.Source2Address.c_str();
	emit signalLog(Address_Source2.c_str());

	std::string Address_Att = "Attenuator address: ";
	Address_Att += heterodyneSettings.Source2Address.c_str();
	emit signalLog(Address_Att.c_str());

	std::string Address_Osc = "Oscilloscope address: ";
	Address_Osc += heterodyneSettings.Source2Address.c_str();
	emit signalLog(Address_Osc.c_str());
	
	emit signalLog("Creating VISA session.");
	error = viOpenDefaultRM(&session);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error in creating VISA: could not locate resources.");
	}
	emit signalLog("VISA session created.");

	emit signalLog("Opening communication channel to PSG1.");
	error = viOpen(session, heterodyneSettings.Source1Address.c_str(), VI_NO_LOCK, 10000, &viPSG1);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error opening channel to PSG1");
	}
	emit signalLog("PSG1 channel opened.");

	emit signalLog("Opening communication channel to PSG2.");
	error = viOpen(session, heterodyneSettings.Source2Address.c_str(), VI_NO_LOCK, 10000, &viPSG2);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error opening channel to PSG2");
	}
	emit signalLog("PSG2 channel opened.");

	emit signalLog("Opening communication channel to Attenuator.");
	error = viOpen(session, heterodyneSettings.AttenuatorAddress.c_str(), VI_NO_LOCK, 10000, &viAttenuator);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error opening channel to Attenuator");
	}
	emit signalLog("Attenuator channel opened.");

	emit signalLog("Opening communication channel to Oscilloscope.");
	error = viOpen(session, heterodyneSettings.OscilloscopeAddress.c_str(), VI_NO_LOCK, 10000, &viOsc);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error opening channel to Oscilloscope");
	}
	emit signalLog("Oscilloscope channel opened.");


	emit signalLog("Asking for PSG1's IDN.");
	error = viPrintf(viPSG1, "*IDN?\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG1.");
	}
	error = viScanf(viPSG1, "%t", buffer);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error reading response from PSG1.");
	}
	std::string source1_idn = "IDN ? -> ";
	source1_idn += buffer;
	emit signalLog(buffer);


	emit signalLog("Asking for PSG2's IDN.");
	error = viPrintf(viPSG2, "*IDN?\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG2.");
	}
	error = viScanf(viPSG2, "%t", buffer);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error reading response from PSG2.");
	}
	std::string source2_idn = "IDN ? -> ";
	source2_idn += buffer;
	emit signalLog(source2_idn.c_str());

	emit signalLog("Asking for Attenuator's IDN.");
	error = viPrintf(viAttenuator, "*IDN?\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to Attenuator.");
	}
	error = viScanf(viAttenuator, "%t", buffer);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error reading response from Attenuator.");
	}
	std::string att_idn = "IDN ? -> ";
	att_idn += buffer;
	emit signalLog(att_idn.c_str());

	emit signalLog("Asking for Oscilloscope's IDN.");
	error = viPrintf(viOsc, "*IDN?\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to Oscilloscope.");
	}
	error = viScanf(viOsc, "%t", buffer);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error reading response from Oscilloscope.");
	}
	std::string osc_idn = "IDN ? -> ";
	osc_idn += buffer;
	emit signalLog(osc_idn.c_str());

	emit signalLog("Resetting Sources.");
	error = viPrintf(viPSG1, ":OUTP 0\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG1.");
	}
	error = viPrintf(viPSG2, ":OUTP 0\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG2.");
	}

	error = viPrintf(viPSG1, ":OUTP:MOD 0\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG1.");
	}
	error = viPrintf(viPSG2, ":OUTP:MOD 0\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG2.");
	}

	error = viPrintf(viPSG1, ":UNIT:POW DBM\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG1.");
	}
	error = viPrintf(viPSG2, ":UNIT:POW DBM\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG2.");
	}

	emit signalLog("Setting sources's output power to defined values.");
	error = viPrintf(viPSG1, ":SOUR:POW:LEV:IMM:AMPL %d\n", heterodyneSettings.source1Amp);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG1.");
	}
	error = viPrintf(viPSG2, ":SOUR:POW:LEV:IMM:AMPL %d\n", heterodyneSettings.source2Amp);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG2.");
	}

	emit signalLog("Setting attenuation to defined value.");
	error = viPrintf(viAttenuator, ":ATT:BANK1:Y %d\n", heterodyneSettings.attenuation);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to Attenuator.");
	}

	emit signalLog("Setting Oscilloscope time range and sample rate.");
	error = viPrintf(viOsc, ":TIM:RANG %f\n", heterodyneSettings.timeRange);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to Oscilloscope.");
	}
	error = viPrintf(viOsc, ":ACQ:SRAT:ANAL %f\n", heterodyneSettings.sampleRate);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to Oscilloscope.");
	}
	error = viPrintf(viOsc, ":RUN\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to Oscilloscope.");
	}

	emit signalLog("Waiting for 5s for Oscilloscope update.");
	sleep(5);

	error = viPrintf(viOsc, ":WAV:POIN?\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to Oscilloscope.");
	}
	error = viScanf(viOsc, "%t", buffer);
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error reading response from Oscilloscope.");
	}

	std::string nPoints_str = buffer;
	int nPoints = std::stoi(nPoints_str);
	emit signalLog("Acquired number of points in Oscilloscope.");
	emit signalLog(nPoints_str.c_str());

	
	char* databuffer = new char[2 * nPoints + 1000];

	double *Y1 = new double[nPoints];
	double *Y2 = new double[nPoints];
	double *Y3 = new double[nPoints];

	std::ofstream datafile(heterodyneSettings.filename);

	datafile << "Frequency" << "," << "I" << "," << "Q" << "\n";

	const double dfreq = (heterodyneSettings.stopFrequency - heterodyneSettings.startFrequency) / (heterodyneSettings.nSteps - 1);
	emit signalRange(heterodyneSettings.stopFrequency, heterodyneSettings.startFrequency);

	const double ifFreq = heterodyneSettings.ifFrequency;


	double I;
	double Q;

	emit signalLog("Starting measurement.");
	error = viPrintf(viPSG1, ":OUTP 1\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG1.");
		goto b11;
	}
	error = viPrintf(viPSG2, ":OUTP 1\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG2.");
		goto b11;
	}

	for (double freq = heterodyneSettings.startFrequency; freq <= heterodyneSettings.stopFrequency; freq += dfreq)
	{

		I = 0;
		Q = 0;

		error = viPrintf(viPSG1, ":FREQ %f\n", freq + ifFreq);
		if (error != VI_SUCCESS)
		{
			emit signalLog("Error sending command to PSG1.");
			goto b11;
		}
		error = viPrintf(viPSG2, ":FREQ %f\n", freq);
		if (error != VI_SUCCESS)
		{
			emit signalLog("Error sending command to PSG2.");
			goto b11;
		}

		usleep(heterodyneSettings.frequency_delay);

		for (int j = 0; j < heterodyneSettings.averages; j++)
		{

			{
				QMutexLocker locker(&m_mutex);
				if (m_stop) 
				{
					goto b11;
				}
			}

			

			error = capture_waveform(viOsc, databuffer, Y1, heterodyneSettings.ChannelI, nPoints);
			if (error != VI_SUCCESS)
			{
				emit signalLog("Error reading Channel I.");
				goto b11;
			}
			error = capture_waveform(viOsc, databuffer, Y2, heterodyneSettings.ChannelSignal, nPoints);
			if (error != VI_SUCCESS)
			{
				emit signalLog("Error reading Channel Signal.");
				goto b11;
			}
			error = capture_waveform(viOsc, databuffer, Y3, heterodyneSettings.ChannelQ, nPoints);
			if (error != VI_SUCCESS)
			{
				emit signalLog("Error reading Channel Q.");
				goto b11;
			}



			double aI = 0;
			double aQ = 0;
			for (int i = 0; i < nPoints; i++)
			{
				aI += Y1[i] * Y2[i] / nPoints;
				aQ += Y3[i] * Y2[i] / nPoints;
			}
			I += aI / heterodyneSettings.averages;
			Q += aQ / heterodyneSettings.averages;

			usleep(heterodyneSettings.acquisition_delay);
		}


		double result = VtodBm(sqrt(I*I + Q * Q));

		emit signalDataPoint(freq, result);

		datafile << freq << "," << I << ',' << Q << "\n";

	}

b11:

	emit signalLog("Finished measurement.");



	error = viPrintf(viPSG1, ":OUTPUT 0\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG1.");
	}
	error = viPrintf(viPSG2, ":OUTPUT 0\n");
	if (error != VI_SUCCESS)
	{
		emit signalLog("Error sending command to PSG2.");
	}


	datafile.close();
b10:
	delete[] Y3;
b9:
	delete[] Y2;
b8:
	delete[] Y1;
b7:
	delete[] databuffer;
b6:
	viClose(viOsc);
b5:
	viClose(viAttenuator);
b4:
	viClose(viPSG2);
b3:
	viClose(viPSG1);
b2:
	viClose(session);
b1:
	delete[] buffer;
}


void HeterodyneThread::run()
{
	if (heterodyneSettings.isSimulation)
	{
		simulate();
	}
	else
	{
		execute();
	}
	
}

void HeterodyneThread::stop()
{
	emit signalLog("Stopping");
	QMutexLocker locker(&m_mutex);
	m_stop = true;
}

HeterodyneThread::HeterodyneThread(std::string dir)
{
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

	const auto acquisition_delay = toml::find<int>(instruments, "Acquisition_delay");
	heterodyneSettings.acquisition_delay = acquisition_delay;

	const auto frequency_change_delay = toml::find<int>(instruments, "Frequency_change_delay");
	heterodyneSettings.frequency_delay = frequency_change_delay;

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