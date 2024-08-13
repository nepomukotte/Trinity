#include "ISiPM.h"

ISiPM::ISiPM(){
	
}
ISiPM::ISiPM(int nPixels){
	gain = std::vector<double>(nPixels, 0);
	optXTalk = std::vector<double>(nPixels, 0);;
	gainTCorrection = std::vector<double>(nPixels, 0);;
	amplToPE = std::vector<double>(nPixels, 0);;
	chgToPE = std::vector<double>(nPixels, 0);;
	
}
ISiPM::~ISiPM(){
	
}

void ISiPM::SetGain(std::vector<double> values){
	gain = values;
}
void ISiPM::SetOptXTalk(std::vector<double> values){
	optXTalk = values;
}
void ISiPM::SetTCorrection(std::vector<double> values){
	gainTCorrection = values;
}
void ISiPM::SetAmplToPE(std::vector<double> values){
	amplToPE = values;
}
void ISiPM::SetChgToPE(std::vector<double> values){
	chgToPE = values;
}

std::vector<double> ISiPM::GetGain(){
	return gain;
}
std::vector<double> ISiPM::GetOptXTalk(){
	return optXTalk;
}
std::vector<double> ISiPM::GetTCorrection(){
	return gainTCorrection;
}
std::vector<double> ISiPM::GetAmplToPE(){
	return amplToPE;
}
std::vector<double> ISiPM::GetChgToPE(){
	return chgToPE;
}