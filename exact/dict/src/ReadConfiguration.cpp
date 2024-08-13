#include "ReadConfiguration.h"

using namespace std;

ReadConfiguration::ReadConfiguration(string cfgFile){
	ifstream *cfgFileStream = new ifstream(cfgFile.c_str());
	string iline;

	pulseExtraction = 0;
	amplitudeExtraction = 0;
	chargeExtraction = 0;
	revolvingTimeWindow = 0;
	isDisplay = 0;
	isFlightMode = 0;


	if(*cfgFileStream){
		while(getline(*cfgFileStream, iline)){
			if (iline.substr(0,1) == "*"){
				ReadLine(iline, cfgFileStream);
			} 
		}
	}
	else{
		cout<<"Read Configuration error: config file not found: "<< cfgFile <<endl;
		exit(-1);
	}

}

ReadConfiguration::~ReadConfiguration(){

}
string ReadConfiguration::BoolToStringAction(Bool_t stage_status){
	string statusText;
	if (stage_status){
		statusText = "YES";
	}
	else{
		statusText = "NO";
	}

	return statusText;
}
void ReadConfiguration::ReadLine(string iline, ifstream *cfgFileStream){
	string i_dump;
	string status;
	istringstream lineStream(iline);
	
	if (iline.find("PULSEEXTRACTION") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>pulseExtraction;

		status = BoolToStringAction(pulseExtraction);

		cout<<"Events will be extracted: "<< status<<endl;
	}
	if (iline.find("EXTRACTAMPLITUDE") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>amplitudeExtraction;

		status = BoolToStringAction(amplitudeExtraction);

		cout<<"Amplitude will be extracted: "<< status<<endl;
	}

	if (iline.find("EXTRACTCHARGE") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>chargeExtraction;

		status = BoolToStringAction(chargeExtraction);

		cout<<"Charge will be extracted: "<< status<<endl;

	}

	if (iline.find("REVOLVINGTIMEWINDOW") < iline.size()){
		lineStream >> i_dump; lineStream >> i_dump;

		lineStream >> revolvingTimeWindow;

		status = BoolToStringAction(revolvingTimeWindow);

		cout << "Will use revolving time window:  "<<status<<endl;
	}

	if (iline.find("SCANTIMEWINDOWSTART") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>timeWindowStart;

		cout<<"Time Window Start: "<< timeWindowStart<<endl;
	}
	if (iline.find("SCANTIMEWINDOWEND") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>timeWindowEnd;

		cout<<"Time Window End: "<< timeWindowEnd<<endl;
	}
	if (iline.find("EXTRACTIONTIMEWINDOW") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>timeWindowExtraction;

		cout<<"Time Window Extraction: "<< timeWindowExtraction<<endl;
	}
	if (iline.find("PEDESTALSAMPLES") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>nSamplesPedestal;

		cout<<"Number of Pedestal Samples: "<< nSamplesPedestal<<endl;
	}
	if (iline.find("NOPIX") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>nPixelsCamera;

		cout<<"Number of Camera Pixels: "<< timeWindowStart<<endl;
	}
	if (iline.find("NOSAMPLES") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>nSamplesADC;

		cout<<"Number of ADC Samples: "<< timeWindowStart<<endl;
	}
	if (iline.find("ADCRANGE") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>nRangeADC;

		cout<<"Number of ADC Samples: "<< nRangeADC<<endl;
	}
	if (iline.find("CALIBRATEPIXELS") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>calibratePixels;
		status = BoolToStringAction(calibratePixels);
		cout<<"Pixels will be calibrated: " << status << endl;
	}
	if (iline.find("CALCULATECALIBRATIONFACTOR") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>calibrationCalculation;

		status = BoolToStringAction(calibrationCalculation);
		cout<<"Pixel calibration calculation: " << status << endl;
	}
	/*if (iline.find("HLEDFILEPREFIX") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>hledFilePrefix;

		cout<<"HLED File for calibration at: "<<hledFilePrefix<<endl;
	}*/
	if (iline.find("ENABLEDISPLAY") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>isDisplay;

		status = BoolToStringAction(isDisplay);
		cout<<"Display: " << status << endl;
	}
	if (iline.find("FLIGHTMODE") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>isFlightMode;

		status = BoolToStringAction(isFlightMode);
		cout<<"Fligt Mode: " << status << endl;
	}
	if (iline.find("ENABLEDISCRIMINATOR") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>isDiscriminator;

		status = BoolToStringAction(isDiscriminator);
		cout<<"Discriminator routine enabled: " << status << endl;
	}
	if (iline.find("ENABLESOFTTHRESHOLD") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>isSoftwareThreshold;

		status = BoolToStringAction(isSoftwareThreshold);
		cout<<"Software Threshold will be used: " << status << endl;
	}
	if (iline.find("SOFTTHRESHOLDLVL") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>softwareThresholdLvl;

		cout<<"Software Threshold level: " << softwareThresholdLvl << endl;
	}
	if (iline.find("GOLDTHRESHOLDLVL") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>goldPlatedThresholdLvl;

		cout<<"Gold Plated Threshold level: " << goldPlatedThresholdLvl << endl;
	}
	if (iline.find("ROI") < iline.size()){
		lineStream>>i_dump; lineStream>>i_dump;
		lineStream>>nRegionOfInterest;

		cout<<"Region of Interest Dimension: " << nRegionOfInterest << endl;
	}
}
