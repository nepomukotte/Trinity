#include "Event.h"

using namespace std;

//ClassImp(Event);
Event::Event(){
	tb_time = 0;
	cobo_time = 0;
	unix_time = 0;
	
	signalValue = vector<vector<Int_t>>(512,vector<Int_t>(512));
	roiMusicID = vector<Int_t>(36);
	roiPixelID = vector<Int_t>(2);
}


Event::~Event(){

}

void Event::SetCoBoTime(ULong64_t rcvTime){
	cobo_time = rcvTime;
}

void Event::SetUNIXTime(ULong64_t rcvTime){
	unix_time = rcvTime;
}

void Event::SetTBTime(ULong64_t rcvTime){
	tb_time = rcvTime;
}

void Event::SetEventType(Int_t rcvEventType){
	eventType = rcvEventType;
}

void Event::SetROIMusicID(vector<Int_t> musicIDs){
	roiMusicID = musicIDs;
}

void Event::SetROIPixelID(vector<Int_t> pixIDs){
	roiPixelID.resize(pixIDs.size());
	roiPixelID = pixIDs;
}

void Event::SetSignalValue(vector<vector<Int_t>> signalTrace){
	signalValue = signalTrace;
	/*for(int i = 0; i< 512; i++){
		cout<<signalValue[0][i]<<endl;
	}*/
}


ULong64_t Event::GetCoBoTime(){
	return cobo_time;
}
ULong64_t Event::GetUNIXTime(){
	return unix_time;
}
ULong64_t Event::GetTBTime(){
	return tb_time;
}
Int_t Event::GetEventType(){
	return eventType;
}
vector<Int_t> Event::GetROIMusicID(){
	return roiMusicID;
}
vector<Int_t> Event::GetROIPixelID(){
	return roiPixelID;
}
vector<vector<Int_t>> Event::GetSignalValue(){
	return signalValue;
}

vector<Int_t> Event::GetSignalValue(int pixID){
	return signalValue[pixID];
}

