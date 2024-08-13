#ifndef PULSE_H
#define PULSE_H
#include <TROOT.h>
#include <ReadConfiguration.h>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

class Pulse {
private:
	Int_t nSamplesPulse;


	Double_t pedestal;
	int amplitude;
	int charge;
	int timePeak;
	int timeFWHM;
	int timeFW;
	float pedestalRMS;
	Double_t timeIntegration;
	vector<Int_t> trace;

	void CalcPedestal(Int_t nSamplesPedestal);
	void CalcAmplitude();
	void CalcCharge(Int_t tWinSt, Int_t tWinEnd, Int_t tWinExt, Bool_t isRevolving);
	void CalcPeakTime();
	void CalcFWHM();
	void CalcFW();
	void FindPeakTime(Int_t tWinSt, Int_t tWinEnd);
	void FindPeak();

public:
	//Constructor
	//Pulse(ReadConfiguration *config, UInt_t *samples);
	Pulse(ReadConfiguration *config, std::vector<Int_t> samples);
	Pulse(vector<Int_t> samples);
	Pulse(vector<int> samples, Int_t tStart, Int_t tEnd, int nSamples = 512, bool isPedestalSub = true);
	Pulse(vector<float> samples, Int_t tStart, Int_t tEnd, int nSamples = 512, bool isPedestalSub = true);
	//Destructor
	~Pulse();



	Double_t GetPedestal();
	float GetPedestalRMS();
	int GetAmplitude();
	int GetCharge();
	int GetTimePeak();
	int GetFWHM();
	int GetFW();		
};
#endif