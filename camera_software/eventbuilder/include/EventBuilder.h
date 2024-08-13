//R__LOAD_LIBRARY(libGET-root.so)
//R__LOAD_LIBRARY(libMultiFrame.so)

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>


#include <mfm/Frame.h>
#include <mfm/Item.h>
#include <mfm/Field.h>
#include <mfm/Serializer.h>
#include <mfm/BitField.h>
#include <mfm/Common.h>
#include <mfm/PrimaryHeader.h>
#include <mfm/Header.h>
#include <mfm/StandardHeader.h>
#include <mfm/BlobHeader.h>
#include <mfm/AbstractFieldContainer.h>
#include <mfm/FrameDictionary.h>

#include <Event.h>


#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include "TBDecoder.cpp"
#include "PixelMap.h"

#define Bifocal     1
#define DiscTest    2
#define Internal    4
#define External    8
#define GPS         16
#define HLED        32

#define ASAD0       0
#define ASAD1       1
#define ASAD2       2
#define ASAD3       3

#define MaxNofChannels			256
#define MaxChannelsPerAget		68
#define ROIMaxNofChannels       36
#define MaxTimeBucket           512
#define HLEDMaxTimeBucket       512
#define BiPedestalTimeBucket    10

using mfm::Frame;

struct STRUCT_Bifocal
{
    ULong64_t tb_time;
    ULong64_t cobo_time;
    ULong64_t unix_time;
    UInt_t ROIMusicID[2];
    UInt_t NofPixels;
    UInt_t ROIPixelID[ROIMaxNofChannels];
    UInt_t ROISignalValue[ROIMaxNofChannels][MaxTimeBucket];
    Float_t Pedestal[MaxNofChannels];
};

struct STRUCT_Background
{
    ULong64_t tb_time;
    ULong64_t cobo_time;
    ULong64_t unix_time;
    UInt_t SignalValue[MaxNofChannels][MaxTimeBucket];
};

struct STRUCT_Hled
{
    ULong64_t tb_time;
    ULong64_t cobo_time;
    ULong64_t unix_time;
    UInt_t SignalValue[MaxNofChannels][HLEDMaxTimeBucket];
};

// Class Definition
class EventBuilder{
private: //member variables

	// pointer to TBDecoder object for accessing triggerboard data
	//TBDecoder *TBEvents;

	// var containing filename
	std::ifstream f_raw;
	std::ifstream f_raw_0;
	std::ifstream f_raw_1;

	// vectors for storing the order of events in a file for each asadID
	// necessary because events are stored sequentially for each asad but not for the file as a whole.
    vector<int> vAsAd0;
    vector<int> vAsAd1;

	// mysterious mfm stuff
	std::auto_ptr<mfm::Frame> curFrame;
	mfm::Item item;
	mfm::Field field;
	mfm::BitField agetIdxField;
	mfm::BitField sampleValueField;

	// current frame header variables from cobo file. Note: 2 frames per event (1 for each asad)
    size_t itemCount_CoBo;
    uint64_t eventTime_CoBo;
    uint32_t eventIdx_CoBo;
    uint8_t asadIdx_CoBo;
    uint32_t eventIdx_CoBo_asad0; //used to confirm data from asad0 and asad1 are paired properly.

    // current event variables from trigger board
    int32_t eventNumber_TB;
    uint64_t eventTime_TB;
    uint8_t triggerSource_TB;
    uint32_t FirstTrigMusic_TB;

	// pointers to TTrees and TFile for storing data
	/*
    TTree *Bifocal_tree;
    TTree *Background_tree;
    TTree *Hled_tree;
    */

    // pixelMapDecoder Variables
    uint16_t sipmID;
    int **pixelMapArray; //stores csv data in heap

	const size_t numChips = 4u;
	const size_t numChannels = 68u;

public:
	// Constructor
	EventBuilder(string formats);
	// Destructor
	~EventBuilder();
private:
	void logTBevent(int eventID, TBDecoder *TBEvents);
	uint16_t FindTrigSource(int eventID);
	int FindPixel(int nx, int ny);
	void FindBin(int iPix,int *nx, int *ny);
	vector<int> FindNeighborPixels(int FirstTrigMusic);
	int FrameCounter(std::string CoBo_filename, int asadIdx);
	void logCoBoEvent(std::auto_ptr<mfm::Frame> &tempFrame);
	/*void Set_Branch_Bifocal(Event *eventBifocal);
	void Set_Branch_Background(Event *eventForced);
	void Set_Branch_Hled(Event *eventHled);*/
	void Write_BifocalEvent(std::string CoBo_filename, Event *eventBifocal, int asadIdx, int event_index, std::vector<std::vector<Int_t>> &signalValue, long startTime);
	void Write_BackgroundEvent(std::string CoBo_filename, Event *eventForced, int asadIdx, int event_index,std::vector<std::vector<Int_t>> &signalValue,long startTime);
	void Write_HledEvent(std::string CoBo_filename, Event *eventHled, int asadIdx, int event_index,std::vector<std::vector<Int_t>> &signalValue,long startTime);
	void Write_TestEvent(std::string CoBo_filename, Event *eventHled, int asadIdx, int event_index,std::vector<std::vector<Int_t>> &signalValue,long startTime);
	//void Write_Tree_to_file();

public:
	// member function definitions
	void mainFlow(std::string CoBo_filename,  std::string filename_tb, std::string outDir);

};
