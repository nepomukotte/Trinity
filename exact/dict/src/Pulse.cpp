#include "Pulse.h"

using namespace std;

/*Pulse::Pulse(ReadConfiguration *config, UInt_t *samples){
	nSamplesPulse = config->nSamplesADC;
	timeIntegration = config->timeWindowExtraction;

	trace = samples;
	amplitude = 0;
	charge = 0;
	timeFWHM = 0;
	timeFW = 0;

	CalcPedestal(config->nSamplesPedestal);
	
	if(config->amplitudeExtraction){
		FindPeakTime(config->timeWindowStart, config->timeWindowEnd);
		FindPeak();
		CalcAmplitude();
	}
	if(config->chargeExtraction){
		CalcCharge(config->timeWindowStart,config->timeWindowEnd,config->timeWindowExtraction,config->revolvingTimeWindow);
	}

}*/

Pulse::Pulse(ReadConfiguration *config, vector<Int_t> samples){
	nSamplesPulse = config->nSamplesADC;
	timeIntegration = config->timeWindowExtraction;

	trace = samples;
	amplitude = 0;
	charge = 0;
	timeFWHM = 0;
	timeFW = 0;
	pedestalRMS = 0;

	CalcPedestal(config->nSamplesPedestal);
	
	if(config->amplitudeExtraction){
		FindPeakTime(config->timeWindowStart, config->timeWindowEnd);
		FindPeak();
		CalcAmplitude();
		//cout<<amplitude<<" ";
	}
	if(config->chargeExtraction){
		CalcCharge(config->timeWindowStart,config->timeWindowEnd,config->timeWindowExtraction,config->revolvingTimeWindow);
	}

}

Pulse::Pulse(vector<Int_t> samples){
	nSamplesPulse = 512;
	timeIntegration = 3;

	trace = samples;
	amplitude = 0;
	charge = 0;
	timeFWHM = 0;
	timeFW = 0;
	pedestalRMS = 0;

	CalcPedestal(100);
	
	FindPeakTime(200, 300);
	FindPeak();
	CalcAmplitude();
		//cout<<amplitude<<" ";
	
	CalcCharge(200,300,3,1);
	
}

Pulse::Pulse(vector<Int_t> samples, Int_t tStart, Int_t tEnd, int nSamples, bool isPedestalSub){
	nSamplesPulse = nSamples;
	timeIntegration = 3;

	trace = samples;
	amplitude = 0;
	charge = 0;
	timeFWHM = 0;
	timeFW = 0;
	pedestalRMS = 0;

	CalcPedestal(100);
	
	FindPeakTime(tStart, tEnd);
	FindPeak();
	if(isPedestalSub){
		CalcAmplitude();
	}
		//cout<<amplitude<<" ";
	
	CalcCharge(tStart,tEnd,3,1);
	
}

Pulse::Pulse(vector<float> samples, Int_t tStart, Int_t tEnd, int nSamples, bool isPedestalSub){
	nSamplesPulse = nSamples;
	timeIntegration = 3;

	trace = std::vector<int>(samples.begin(), samples.end());
	amplitude = 0;
	charge = 0;
	timeFWHM = 0;
	timeFW = 0;
	pedestalRMS = 0;

	CalcPedestal(100);
	
	FindPeakTime(tStart, tEnd);
	FindPeak();
	if(isPedestalSub){
		CalcAmplitude();
	}
	
		//cout<<amplitude<<" ";
	
	CalcCharge(tStart,tEnd,3,1);
	
}






Pulse::~Pulse(){

}

void Pulse::FindPeakTime(Int_t tWinSt, Int_t tWinEnd){
	Double_t aCalc, aTmp;
	Double_t pkTmp = tWinSt;

	aCalc = 5000;

	
	for(int i = tWinSt; i<tWinEnd; i++){
		aTmp = trace[i];
		if (aTmp<aCalc){
			aCalc = aTmp;
			pkTmp = i;
		}
	}

	timePeak = pkTmp;
}

//Calculate desired properties of Pulse
void Pulse::FindPeak(){
	amplitude = trace[(int)timePeak];
}

void Pulse::CalcPedestal(Int_t nSamplesPedestal){
	Double_t pCalc = 0;

	for(int i = 0; i<nSamplesPedestal; i++){
		pCalc += trace[i];
	}

	pCalc = pCalc/nSamplesPedestal;

	pedestal = pCalc;

	for(int i = 0; i<nSamplesPedestal; i++){
		pedestalRMS += (trace[i] - pedestal)*(trace[i] - pedestal);
	}
	pedestalRMS = sqrt(pedestalRMS/(float)nSamplesPedestal);
	//cout<<pedestalRMS<<endl;
}


void Pulse::CalcAmplitude(){
	amplitude = pedestal - amplitude;
	if(amplitude<0){
		amplitude = 0;
	}
//=======
//	amplitude = pedestal - amplitude;
//	if(amplitude<0){
//		amplitude = 0;
//	}
//>>>>>>> main
}

void Pulse::CalcCharge(Int_t tWinSt, Int_t tWinEnd, Int_t tWinExt, Bool_t isRevolving){
	if(isRevolving){
		FindPeakTime(tWinSt, tWinEnd);
		for(int i = timePeak - tWinExt/2; i<timePeak+tWinExt/2+tWinExt%2; i++){
			charge += (pedestal - trace[i]);
		}

	}
	else{
		for(int i = tWinSt; i<tWinSt+tWinExt; i++){
			charge += (pedestal - trace[i]);
		}
	}

	if(charge<0){
		charge=0;
	}

}

void Pulse::CalcFWHM(){
	
}

void Pulse::CalcFW(){

}

//Getters
Double_t Pulse::GetPedestal(){
	return pedestal;
}

float Pulse::GetPedestalRMS(){
	return pedestalRMS;
}

int Pulse::GetAmplitude(){
	return amplitude;
}

int Pulse::GetCharge(){
	return charge;
}

int Pulse::GetTimePeak(){
	return timePeak;
}

int Pulse::GetFWHM(){
	return timeFWHM;
}

int Pulse::GetFW(){
	return timeFW;
}
