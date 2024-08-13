#ifndef EXTRACTEDDATA_H
#define EXTRACTEDDATA_H

#include <TROOT.h>
#include "ReadConfiguration.h"

using namespace std;

class ExtractedData{

public:
	ExtractedData();
	~ExtractedData();

	//Setters
	void SetAmplitude(unsigned short ampl);
	void SetCharge(unsigned short chg);
	void SetTimePeak(unsigned short tmPk);
	//void SetTimeFWHM(unsigned short tmFWHM);
	//void SetTimeFW(unsigned short tmFW);
	//void SetCorrFactor(float crFctr);
	void SetPedestal(unsigned short rcvPedestal);
	void SetPedestalRMS(float rcvPedestalRMS);
	//void SetEventTime(float rcvEventTime);

	unsigned short GetAmplitude();
	unsigned short GetCharge();
	unsigned short GetTimePeak();
	unsigned short GetPedestal();
	float GetPedestalRMS();

	void Copy(ExtractedData *ex);
	//float GetEventTime();
	//unsigned short GetTimeFWHM();
	//unsigned short GetTimeFW();
	//float GetCorrFactor();



protected:
	unsigned short amplitude;
	unsigned short charge;
	unsigned short timePeak;
	unsigned short pedestal;
	float pedestalRMS;
	//float eventTime;
	//unsigned short timeFWHM;
	//unsigned short timeFW;

	//float corrFactor;

//ClassDef(ExtractedData,1);
};
#endif