//! TriggerRead reading a configuration file for the trigger analysis
#ifndef READCONFIG_H
#define READCONFIG_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <TROOT.h>

using namespace std;

class ReadConfiguration{
public:
	ReadConfiguration(string cfgFile);
	~ReadConfiguration();

	Bool_t pulseExtraction;
	Bool_t amplitudeExtraction;
	Bool_t chargeExtraction;
	Bool_t revolvingTimeWindow;
	Bool_t calibrationCalculation;
	Bool_t calibratePixels;
	Bool_t isDisplay;
	Bool_t isFlightMode;
	Bool_t isDiscriminator;
	Bool_t isSoftwareThreshold;

	Int_t timeWindowStart;
	Int_t timeWindowEnd;
	Int_t timeWindowExtraction;

	Int_t nSamplesPedestal;
	Int_t nPixelsCamera;
	Int_t nSamplesADC;
	Int_t nRangeADC;
	Int_t nRegionOfInterest;

	Int_t softwareThresholdLvl;
	Int_t goldPlatedThresholdLvl;

protected:
	void ReadLine(string iline, ifstream *inFileStream);
	string BoolToStringAction(Bool_t stage_status);
};

#endif