#include "ExtractedData.h"

using namespace std;

//ClassImp(ExtractedData);
ExtractedData::ExtractedData(){
	amplitude = 0;
	charge = 0;
	timePeak = 0;
	pedestal = 0;
	pedestalRMS = 0;
	//timeFWHM = 0;
	//timeFW = 0;
	//corrFactor = 1;
}

ExtractedData::~ExtractedData(){

}
//Setters
void ExtractedData::SetAmplitude(unsigned short ampl){
	amplitude = ampl;
}

void ExtractedData::SetCharge(unsigned short chg){
	charge = chg;
}

void ExtractedData::SetTimePeak(unsigned short tmPk){
	timePeak = tmPk;
}

void ExtractedData::SetPedestalRMS(float pdrms){
	pedestalRMS = pdrms;
}

/*void ExtractedData::SetTimeFWHM(int tmFWHM){
	timeFWHM = tmFWHM;
}

void ExtractedData::SetTimeFW(int tmFW){
	timeFW = tmFW;
}*/

/*
void ExtractedData::SetCorrFactor(float crFctr){
	corrFactor = crFctr;
}
*/
/*void ExtractedData::SetEventTime(float rcvEventTime){
	eventTime = rcvEventTime;
}*/

void ExtractedData::SetPedestal(unsigned short rcvPedestal){
	pedestal = rcvPedestal;
}

//Getters

unsigned short ExtractedData::GetAmplitude(){
	return amplitude;
}
unsigned short ExtractedData::GetCharge(){
	return charge;
}
unsigned short ExtractedData::GetTimePeak(){
	return timePeak;
}
/*unsigned short ExtractedData::GetTimeFWHM(){
	return timeFWHM;
}
unsigned short ExtractedData::GetTimeFW(){
	return timeFW;
}*/
/*
float ExtractedData::GetCorrFactor(){
	return corrFactor;
}
*/
/*float ExtractedData::GetEventTime(){
	return eventTime;
}*/

unsigned short ExtractedData::GetPedestal(){
	return pedestal;
}

float ExtractedData::GetPedestalRMS(){
	return pedestalRMS;
}

void ExtractedData::Copy(ExtractedData *ex){
	amplitude = ex->GetAmplitude();
	charge = ex->GetCharge();
	pedestal = ex->GetPedestal();
	timePeak = ex->GetTimePeak();
	pedestalRMS = ex->GetPedestalRMS();
}
