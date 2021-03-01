#pragma once

struct MeasurementSetting {
	double startFrequency;
	double stopFrequency;
	int nSteps;
	int attenuation;
	int averages;
	double timeRange;
	double sampleRate;
	double ifFrequency;
	std::string filename;
	std::string Source1Address;
	std::string Source2Address;
	std::string AttenuatorAddress;
	std::string OscilloscopeAddress;
	int acquisition_delay;
	int frequency_delay;
	int source1Amp;
	int source2Amp;
	bool isSimulation;
	int ChannelI;
	int ChannelQ;
	int ChannelSignal;
};