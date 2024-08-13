#ifndef EVENT_H 
#define EVENT_H

#include <TROOT.h>
#include "ReadConfiguration.h"
#include <vector>
using namespace std;

class Event {
public:
	Event();
	~Event();

	void SetCoBoTime(ULong64_t rcvTime);
	void SetUNIXTime(ULong64_t rcvTime);
	void SetTBTime(ULong64_t rcvTime);
	void SetEventType(Int_t rcvEventType);
	void SetROIMusicID(vector<Int_t> musicIDs);
	void SetROIPixelID(vector<Int_t> pixIDs);
	void SetSignalValue(vector<vector<Int_t>> signalTrace);

	ULong64_t GetCoBoTime();
	ULong64_t GetUNIXTime();
	ULong64_t GetTBTime();
	Int_t GetEventType();
	vector<Int_t> GetROIMusicID();
	vector<Int_t> GetROIPixelID();
	vector<vector<Int_t>> GetSignalValue();
	vector<Int_t> GetSignalValue(int pixID);
	

protected:
	ULong64_t tb_time;
	ULong64_t cobo_time;
	ULong64_t unix_time;
	Int_t eventType;

	vector<vector<Int_t>> signalValue;
	vector<Int_t> roiMusicID;
	vector<Int_t> roiPixelID;

//ClassDef(Event,1);
};
#endif