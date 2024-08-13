#include "GoldPlated.h"

using namespace std;
GoldPlated::GoldPlated(){
	signalValue = vector<unsigned short>();
}

GoldPlated::GoldPlated(int traceSize){
	signalValue = vector<unsigned short>(traceSize);
}

GoldPlated::~GoldPlated(){

}

void GoldPlated::SetTrace(vector<int> trace, int tStart){
	signalValue = vector<unsigned short>(trace.begin()+tStart,trace.begin()+tStart+signalValue.size());
	
}

void GoldPlated::CopyFromExtractedData(ExtractedData *ex){
	//pxIDs = bf->GetROIPixelIDs();
	//musicIDs = bf->GetTrigMUSICIDs();

	amplitude = ex->GetAmplitude();
	charge = ex->GetCharge();
	timePeak = ex->GetTimePeak();
	pedestal = ex->GetPedestal();
	pedestalRMS = ex->GetPedestalRMS();
	//eventTime = ex->GetEventTime();

}

vector<unsigned short> GoldPlated::GetTrace(){
	vector<unsigned short> pxTrace = signalValue;
	return pxTrace;
}

// vector<vector<unsigned short>> GoldPlated::GetAllTraces(){
// 	return signalValue;
// }