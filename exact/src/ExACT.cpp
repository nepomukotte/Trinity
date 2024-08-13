#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TList.h>
#include <TKey.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TTimer.h>
#include <TMath.h>
#include <exception>
#include <fstream>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TDatime.h>
#include <chrono>
#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>

// libExACT objects
#include <ExtractedData.h>
#include <ReadConfiguration.h>
#include <Event.h>
#include <GoldPlated.h>
#include <BiFocal.h>
#include <Pulse.h>

// ExACT Objects
//#include "Pulse.h"
#include <IPlotTools.h>



using namespace std;

bool CheckTriggerGeometry(unsigned short timePeak_1, unsigned short timePeak_2);
bool CheckIsAboveThreshold(unsigned short amplitude_1, unsigned short amplitude_2, int threshold);

bool CheckFileExists(string fileName){
	bool isFileExists;
	ifstream f(fileName.c_str());
	isFileExists = f.good();

	return isFileExists;
}

void StarlinkFileExtraction(ReadConfiguration *readConf,string dataFilePrefix){
	string dataFileName = dataFilePrefix + ".root";
	
	string extractedFileName = dataFilePrefix + "_Starlink"+".root";
	//string treeName;
	Int_t nEventsHLED, nEventsForced, nEventsBiFocal, nEventsTest;

	nEventsHLED = 0;
	nEventsForced = 0;
	nEventsBiFocal = 0;
	nEventsTest = 0;

	TFile *fileIn = new TFile(dataFileName.c_str(),"READ");
	TTree *treeInHLED = (TTree*)fileIn->Get("HLED");
	TTree *treeInForced = (TTree*)fileIn->Get("Forced");
	TTree *treeInBiFocal = (TTree*)fileIn->Get("BiFocal");
	TTree *treeInTest = (TTree*)fileIn->Get("Test");


	Event *eventBifocal = new Event();
	Event *eventTest = new Event();

	//treeInHLED->SetBranchAddress("Events", &eventHLED);
	//treeInForced->SetBranchAddress("Events", &eventForced);
	treeInBiFocal->SetBranchAddress("Events", &eventBifocal);
	treeInTest->SetBranchAddress("Events",&eventTest);

	// TTree *treeH= (TTree*)fileIn->Get(treeName.c_str());



	TFile *fileOut = new TFile	(extractedFileName.c_str(),"RECREATE");
	TTree *treeOutHLED = 0;
	treeOutHLED = treeInHLED->CloneTree(5);
	TTree *treeOutForced = 0;
	treeOutForced = treeInForced->CloneTree();

	TTree *treeOutBifocal = treeInBiFocal->CloneTree(0);
	TTree *treeOutTest = treeInTest->CloneTree(0);
	//TTree *treeOutCorrectionFactor = new TTree("Corrections", "Calibration Factors");

	//vector<float> corrFactor = vector<float>(readConf->nPixelsCamera);



	
	nEventsBiFocal = treeInBiFocal->GetEntries();
	nEventsTest = treeInTest->GetEntries();

	//cout<<"Here"<<endl;

	Pulse *pulse;

	Pulse *pulseLeft, *pulseRight;

	vector<Int_t> musicIDsROI = vector<Int_t>(2);
	Int_t maxChg = 0;
	int maxPxID;
	// Extract all BiFocal Events in file
	for(int i = 0; i<nEventsBiFocal; i++){
		fileIn->cd();
		treeInBiFocal->GetEntry(i);
		musicIDsROI = eventBifocal->GetROIMusicID();
		//cout<<"Now here"<<endl;
		for(int j = musicIDsROI[0]*8; j<musicIDsROI[0]*8+8; j++){
			pulse = new Pulse(readConf,eventBifocal->GetSignalValue(j));
			if(pulse->GetAmplitude()>maxChg){
				maxPxID = j;
				maxChg = pulse->GetAmplitude();
			}

		}

		//cout<<"Out of the woods"<<endl;
		pulseLeft = new Pulse(readConf,eventBifocal->GetSignalValue(maxPxID));
		pulseRight = new Pulse(readConf,eventBifocal->GetSignalValue(maxPxID+8));
		if(CheckTriggerGeometry(pulseLeft->GetTimePeak(),pulseRight->GetTimePeak())){
			if(CheckIsAboveThreshold(pulseLeft->GetAmplitude(), 
				pulseRight->GetAmplitude(),
				readConf->softwareThresholdLvl)){
					fileOut->cd();
					treeOutBifocal->Fill();
					//cout<<"Event: "<<i<<" survived"<<endl;
			}
		}
		//cout<<endl;
		if (i%50 ==0){
			cout << "Extracting Event: "<<i<<" from BiFocal tree"<<endl;
		}

	}

	maxChg = 0;

	for(int i = 0; i<nEventsTest; i++){
		fileIn->cd();
		treeInTest->GetEntry(i);
		musicIDsROI = eventTest->GetROIMusicID();
		for(int j = musicIDsROI[0]*8; j<musicIDsROI[0]*8+8; j++){
			pulse = new Pulse(readConf,eventTest->GetSignalValue(j));
			if(pulse->GetAmplitude()>maxChg){
				maxPxID = j;
				maxChg = pulse->GetAmplitude();
			}
			pulseLeft = new Pulse(readConf,eventTest->GetSignalValue(maxPxID));
			pulseRight = new Pulse(readConf,eventTest->GetSignalValue(maxPxID+8));
			if(CheckIsAboveThreshold(pulseLeft->GetAmplitude(),
				pulseRight->GetAmplitude(),
				readConf->softwareThresholdLvl)){
					fileOut->cd();
					treeOutTest->Fill();
			}
		}
		if (i%50 ==0){
			cout << "Extracting Event: "<<i<<" from Test tree"<<endl;
		}

	}

	fileOut->cd();

	treeOutHLED->Fill();
	treeOutTest->Fill();

	fileOut->Write();
	fileOut->Close();
	fileIn->Close();
	//treeOut->Write();
	//delete pulse;
	//delete pulseLeft;
	//delete pulseRight;
	cout<<"Done"<<endl;
}

void PulseExtraction(ReadConfiguration *readConf,string dataFilePrefix){
	string dataFileName = dataFilePrefix + ".root";
	
	string extractedFileName = dataFilePrefix + "_Extracted"+".root";
	//string treeName;
	Int_t nEventsHLED, nEventsForced, nEventsBiFocal, nEventsTest;

	nEventsHLED = 0;
	nEventsForced = 0;
	nEventsBiFocal = 0;
	nEventsTest = 0;

	TFile *fileIn = new TFile(dataFileName.c_str(),"READ");
	TTree *treeInHLED = (TTree*)fileIn->Get("HLED");
	TTree *treeInForced = (TTree*)fileIn->Get("Forced");
	TTree *treeInBiFocal = (TTree*)fileIn->Get("BiFocal");
	TTree *treeInTest = (TTree*)fileIn->Get("Test");

	Event *eventHLED = new Event();
	Event *eventForced = new Event();
	Event *eventBifocal = new Event();
	Event *eventTest = new Event();

	treeInHLED->SetBranchAddress("Events", &eventHLED);
	treeInForced->SetBranchAddress("Events", &eventForced);
	treeInBiFocal->SetBranchAddress("Events", &eventBifocal);
	treeInTest->SetBranchAddress("Events",&eventTest);

	// TTree *treeH= (TTree*)fileIn->Get(treeName.c_str());



	TFile *fileOut = new TFile	(extractedFileName.c_str(),"RECREATE");
	TTree *treeOutHLED = new TTree("HLED","Extracted Data");
	TTree *treeOutForced = new TTree("Forced","Extracted Data");
	TTree *treeOutBifocal = new TTree("BiFocal","Extracted Data");
	TTree *treeOutTest = new TTree("Test","Extracted Data");
	//TTree *treeOutCorrectionFactor = new TTree("Corrections", "Calibration Factors");

	//vector<float> corrFactor = vector<float>(readConf->nPixelsCamera);

	ExtractedData *extractedDataHLED[readConf->nPixelsCamera];
	ExtractedData *extractedDataForced[readConf->nPixelsCamera];
	ExtractedData *extractedDataBifocal[readConf->nPixelsCamera];
	ExtractedData *extractedDataTest[readConf->nPixelsCamera];

	unsigned long long triggerTimeHLED, triggerTimeForced, triggerTimeBifocal, triggerTimeTest;
	triggerTimeHLED = 0;
	triggerTimeForced = 0;
	triggerTimeBifocal = 0;
	triggerTimeTest = 0;

	BiFocal *bifocalInfo = 0;
	bifocalInfo = new BiFocal();

	BiFocal *bifocalInfoTest = 0;
	bifocalInfoTest = new BiFocal();

	for (int i = 0; i<readConf->nPixelsCamera; i++){
		extractedDataHLED[i] = 0;
		extractedDataForced[i] = 0;
		extractedDataBifocal[i] = 0;
		extractedDataTest[i] = 0;

		extractedDataHLED[i] = new ExtractedData();
		extractedDataForced[i] = new ExtractedData();
		extractedDataBifocal[i] = new ExtractedData();
		extractedDataTest[i] = new ExtractedData();	
	}

	for (int i = 0; i<readConf->nPixelsCamera; i++){
		treeOutHLED->Branch(TString::Format("Pixel[%d]",i),
			"ExtractedData",&extractedDataHLED[i],64000,0);
		treeOutForced->Branch(TString::Format("Pixel[%d]",i),
			"ExtractedData",&extractedDataForced[i],64000,0);
		treeOutBifocal->Branch(TString::Format("Pixel[%d]",i),
			"ExtractedData",&extractedDataBifocal[i],64000,0);
		treeOutTest->Branch(TString::Format("Pixel[%d]",i),
			"ExtractedData",&extractedDataTest[i],64000,0);
		//treeOutCorrectionFactor->Branch(TString::Format("Pixel[%d]",i),&corrFactor[i]);
		//cout<<&extractedData[i]<<endl;
	}
	treeOutHLED->Branch("TriggerTime", &triggerTimeHLED);
	treeOutForced->Branch("TriggerTime", &triggerTimeForced);
	treeOutBifocal->Branch("TriggerTime", &triggerTimeBifocal);
	treeOutTest->Branch("TriggerTime", &triggerTimeTest);

	treeOutBifocal->Branch("TriggerRegion","BiFocal",bifocalInfo,64000,0);
	treeOutTest->Branch("TriggerRegion","BiFocal",bifocalInfoTest,64000,0);


	nEventsHLED = treeInHLED->GetEntries();
	nEventsForced = treeInForced->GetEntries();
	nEventsBiFocal = treeInBiFocal->GetEntries();
	nEventsTest = treeInTest->GetEntries();


	Pulse *pulse;

	// Extract all of the HLED events in the included file
	for(int i = 0; i<nEventsHLED; i++){
		fileIn->cd();
		treeInHLED->GetEntry(i);
		for(int j = 0; j<readConf->nPixelsCamera; j++){
			pulse = new Pulse(readConf,eventHLED->GetSignalValue(j));
			if(readConf->amplitudeExtraction){
				extractedDataHLED[j]->SetAmplitude(pulse->GetAmplitude());
				extractedDataHLED[j]->SetTimePeak(pulse->GetTimePeak());
			}
			if(readConf->chargeExtraction){
				extractedDataHLED[j]->SetCharge(pulse->GetCharge());
			}
			if(j==322){
				//cout<<"Charge"<<pulse->GetAmplitude()<<endl;;
			}
			extractedDataHLED[j]->SetPedestal(round(pulse->GetPedestal()));
			extractedDataHLED[j]->SetPedestalRMS((pulse->GetPedestalRMS()));
			//extractedDataHLED[j]->SetEventTime(eventHLED->GetTBTime());
			/*if(readConf->revolvingTimeWindow){
				extractedData[j]->SetTimePeak(pulse->GetTimePeak());
			}*/
		}

		triggerTimeHLED = eventHLED->GetTBTime();

		if (i%50 ==0){
			cout << "Extracting Event: "<<i<<" from HLED tree"<<endl;
		}
		fileOut->cd();
		treeOutHLED->Fill();

	}

	// Extract all Forced trigger events in file
	for(int i = 0; i<nEventsForced; i++){
		fileIn->cd();
		treeInForced->GetEntry(i);
		for(int j = 0; j<readConf->nPixelsCamera; j++){
			pulse = new Pulse(readConf,eventForced->GetSignalValue(j));
			if(readConf->amplitudeExtraction){
				extractedDataForced[j]->SetAmplitude(pulse->GetAmplitude());
				extractedDataForced[j]->SetTimePeak(pulse->GetTimePeak());
			}
			if(readConf->chargeExtraction){
				extractedDataForced[j]->SetCharge(pulse->GetCharge());
			}
			extractedDataForced[j]->SetPedestal(round(pulse->GetPedestal()));
			extractedDataForced[j]->SetPedestalRMS((pulse->GetPedestalRMS()));
			//extractedDataForced[j]->SetEventTime(eventHLED->GetTBTime());
			/*if(readConf->revolvingTimeWindow){
				extractedData[j]->SetTimePeak(pulse->GetTimePeak());
			}*/
		}
		triggerTimeForced = eventForced->GetTBTime();
		if (i%50 ==0){
			cout << "Extracting Event: "<<i<<" from Forced tree"<<endl;
		}
		fileOut->cd();
		treeOutForced->Fill();

	}

	// Extract all BiFocal Events in file
	for(int i = 0; i<nEventsBiFocal; i++){
		fileIn->cd();
		treeInBiFocal->GetEntry(i);
		for(int j = 0; j<readConf->nPixelsCamera; j++){
			pulse = new Pulse(readConf,eventBifocal->GetSignalValue(j));
			//sleep(3);
			//cout<<pulse->GetAmplitude()<<" "<<endl;

			if(readConf->amplitudeExtraction){
				
				extractedDataBifocal[j]->SetAmplitude(pulse->GetAmplitude());
				extractedDataBifocal[j]->SetTimePeak(pulse->GetTimePeak());
			}
			if(readConf->chargeExtraction){
				extractedDataBifocal[j]->SetCharge(pulse->GetCharge());
			}
			extractedDataBifocal[j]->SetPedestal(round(pulse->GetPedestal()));
			extractedDataBifocal[j]->SetPedestalRMS((pulse->GetPedestalRMS()));
			//extractedDataBifocal[j]->SetEventTime(eventHLED->GetTBTime());
			
		}

		//cout<<endl;
		
		triggerTimeBifocal = eventBifocal->GetTBTime();
		bifocalInfo->SetROIPixelIDs(eventBifocal->GetROIPixelID());
		bifocalInfo->SetTrigMUSICIDs(eventBifocal->GetROIMusicID());
		if (i%50 ==0){
			cout << "Extracting Event: "<<i<<" from BiFocal tree"<<endl;
		}
		fileOut->cd();
		treeOutBifocal->Fill();

	}

	for(int i = 0; i<nEventsTest; i++){
		fileIn->cd();
		treeInTest->GetEntry(i);
		for(int j = 0; j<readConf->nPixelsCamera; j++){
			pulse = new Pulse(readConf,eventTest->GetSignalValue(j));
			//sleep(3);
			//cout<<pulse->GetAmplitude()<<" "<<endl;

			if(readConf->amplitudeExtraction){
				
				extractedDataTest[j]->SetAmplitude(pulse->GetAmplitude());
				extractedDataTest[j]->SetTimePeak(pulse->GetTimePeak());
			}
			if(readConf->chargeExtraction){
				extractedDataTest[j]->SetCharge(pulse->GetCharge());
			}
			extractedDataTest[j]->SetPedestal(round(pulse->GetPedestal()));
			extractedDataTest[j]->SetPedestalRMS((pulse->GetPedestalRMS()));
			//extractedDataBifocal[j]->SetEventTime(eventHLED->GetTBTime());
			
		}

		//cout<<endl;
		
		triggerTimeTest = eventTest->GetTBTime();
		bifocalInfoTest->SetROIPixelIDs(eventTest->GetROIPixelID());
		bifocalInfoTest->SetTrigMUSICIDs(eventTest->GetROIMusicID());
		if (i%50 ==0){
			cout << "Extracting Event: "<<i<<" from Test tree"<<endl;
		}
		fileOut->cd();
		treeOutTest->Fill();

	}
	
	fileOut->Write();
	fileOut->Close();
	fileIn->Close();
	//treeOut->Write();
	cout<<"Done"<<endl;
}

void LoadCorrectionFactors(ReadConfiguration *readConf, string dataFilePrefix, vector<float> &correctionFactor){
	string correctionFactorFileName = dataFilePrefix + "_CorrectionFactor.txt";
	
	string line;
	ifstream fileCorrectionFactor(correctionFactorFileName.c_str());
	int i = 0;
	while(getline(fileCorrectionFactor,line)){
		//cout<<"Line: "<<line<<endl;
		correctionFactor[i] = stof(line);
		//cout<<"Line: "<<correctionFactor[i]<<endl;
		i++;
	}

	fileCorrectionFactor.close();
}

bool CheckTriggerGeometry(unsigned short timePeak_1, unsigned short timePeak_2){
	bool isGeometryCorrect = false;

	int tolerance = 5;
	//cout<<"Time: "<<timePeak_1<<" "<<timePeak_2<<endl;
	if(TMath::Abs(timePeak_1-timePeak_2)<=tolerance){
		isGeometryCorrect = true;
	}else{
		isGeometryCorrect = false;
	}

	return isGeometryCorrect;
}

bool CheckIsAboveThreshold(unsigned short amplitude_1, unsigned short amplitude_2, int threshold){
	bool isAboveThreshold = false;
	int ADCPerPe = 1;
	//cout<<amplitude_1<<" "<<amplitude_2<<" "<<threshold<<endl;
	if (amplitude_1 >= threshold*ADCPerPe && amplitude_2 >= threshold*ADCPerPe){
		isAboveThreshold = true;
	}
	else{
		isAboveThreshold = false;
	}

	return isAboveThreshold;
}

void LoadTracesForGoldPlated(ReadConfiguration *readConf, int eventIdx, GoldPlated **gP, string dataFilePrefix){
	string dataFileName = dataFilePrefix + ".root";

	TFile *fileEvents = new TFile(dataFileName.c_str(),"READ");

	TTree *treeEventsBifocal = (TTree*)fileEvents->Get("BiFocal");

	Event *ev = new Event();

	treeEventsBifocal->SetBranchAddress("Events",&ev);

	treeEventsBifocal->GetEntry(eventIdx);
	for(int j = 0; j<readConf->nPixelsCamera; j++){
		gP[j]->SetTrace(ev->GetSignalValue(j), readConf->timeWindowStart);
	}
	fileEvents->Close();
}

void CopyDataFromExtractedBiFocal(GoldPlated *gp, ExtractedData *bf){
	gp->CopyFromExtractedData(bf);
	
}
void CopyDataFromExtractedBiFocal(ExtractedData *ex, ExtractedData *bf){
	ex->Copy(bf);
	
}

void EventDiscrimination(ReadConfiguration *readConf, string dataFilePrefix){
	string extractedFileName = dataFilePrefix + "_Extracted.root";
	string downloadFileName = dataFilePrefix + "_Download.root";
	
	//cout<<"Here"<<endl;

	TFile *fileExtractedData = new TFile(extractedFileName.c_str(), "READ");

	TTree *treeInHLED = (TTree*)fileExtractedData->Get("HLED");
	TTree *treeInBiFocal = (TTree*)fileExtractedData->Get("BiFocal");
	TTree *treeInForced = (TTree*)fileExtractedData->Get("Forced");
	TTree *treeInTest = (TTree*)fileExtractedData->Get("Test");

	Int_t nEventsBiFocal = 0;
	Int_t nEventsTest = 0;

	nEventsBiFocal = treeInBiFocal->GetEntries();
	nEventsTest = treeInTest->GetEntries();

	ExtractedData *extractedDataBifocal[readConf->nPixelsCamera];
	ExtractedData *extractedDataTest[readConf->nPixelsCamera];
	BiFocal *bifocalInfo = 0;
	BiFocal *bifocalInfoTest = 0;
	unsigned long long triggerTimeBifocal = 0;
	unsigned long long triggerTimeTest = 0;

	for (int i = 0; i<readConf->nPixelsCamera; i++){
		extractedDataBifocal[i] = 0;
		extractedDataBifocal[i] = new ExtractedData();
		treeInBiFocal->SetBranchAddress(TString::Format("Pixel[%d]",i),&extractedDataBifocal[i]);
		
		extractedDataTest[i] = 0;
		extractedDataTest[i] = new ExtractedData();
		treeInTest->SetBranchAddress(TString::Format("Pixel[%d]",i),&extractedDataTest[i]);
	}

	bifocalInfo = new BiFocal();
	bifocalInfoTest = new BiFocal();

	treeInTest->SetBranchAddress("TriggerRegion",&bifocalInfoTest);
	treeInTest->SetBranchAddress("TriggerTime",&triggerTimeTest);

	treeInBiFocal->SetBranchAddress("TriggerRegion",&bifocalInfo);
	treeInBiFocal->SetBranchAddress("TriggerTime", &triggerTimeBifocal);


	


	//cout<<"Here Again"<<endl;
	TFile *fileDownloadReady = new TFile(downloadFileName.c_str(),"RECREATE");

	vector<unsigned short> pixIDsROI;
	vector<unsigned short> musicIDsROI;
	
	//vector<float*> memCorrFactor = vector<float*>(readConf->nPixelsCamera); 
	


	TTree *treeOutHLED = 0;//new TTree("HLED","Extracted Data");
	TTree *treeOutBiFocal = new TTree("BiFocal","Extracted Data");
	TTree *treeOutGoldPlated = new TTree("GoldPlated","Extracted Data");
	TTree *treeOutGoldPlatedTest = new TTree("GoldPlatedTest","Extracted Data");
	TTree *treeOutForced = 0;//new TTree("Forced","Extracted Data");
	TTree *treeOutTest = new TTree("Test","ExtractedData");
	TTree *treeOutCorrectionFactor = new TTree("CorrectionFactor","Correction Factors");

	vector<float> *correctionFactorSave = 0;
	vector<float> correctionFactors = vector<float>(readConf->nPixelsCamera,1.0);

	ExtractedData *extractedDataBifocalOut[readConf->nPixelsCamera];
	ExtractedData *extractedDataTestOut[readConf->nPixelsCamera];
	GoldPlated *goldPlatedOut[readConf->nPixelsCamera];
	GoldPlated *goldPlatedOutTest[readConf->nPixelsCamera];

	BiFocal *bifocalInfoBifocal = 0;
	BiFocal *bifocalInfoTestOut =0;
	bifocalInfoTestOut = new BiFocal();
	bifocalInfoBifocal = new BiFocal();


	BiFocal *bifocalInfoGoldPlated = 0;
	BiFocal *bifocalInfoGoldPlatedTest = 0;
	bifocalInfoGoldPlatedTest = new BiFocal();
	bifocalInfoGoldPlated = new BiFocal();

	unsigned long long  triggerTimeBifocalOut, triggerTimeGoldPlatedOut, triggerTimeTestOut, triggerTimeGoldPlatedTestOut;

	triggerTimeBifocalOut = 0;
	triggerTimeGoldPlatedOut = 0;

	triggerTimeTestOut = 0;
	triggerTimeGoldPlatedTestOut = 0;

	TString branchNameE, branchNameG, branchNameC, branchNameET, branchNameGT;
	//cout << "Here 4" <<endl;
	for(int i = 0; i<readConf->nPixelsCamera; i++){
		extractedDataBifocalOut[i] = 0;
		extractedDataTestOut[i] = 0;
		goldPlatedOut[i] = 0;
		goldPlatedOutTest[i] = 0;

		extractedDataBifocalOut[i] = new ExtractedData();
		extractedDataTestOut[i] = new ExtractedData();

		goldPlatedOutTest[i] = new GoldPlated(readConf->timeWindowEnd-readConf->timeWindowStart);
		goldPlatedOut[i] = new GoldPlated(readConf->timeWindowEnd-readConf->timeWindowStart);
	}
	for(int i = 0; i<readConf->nPixelsCamera; i++){
		branchNameE = TString::Format("Pix[%d]",i);
		treeOutBiFocal->Branch(branchNameE,
			"ExtractedData",&extractedDataBifocalOut[i],64000,0);
	}
	for(int i = 0; i<readConf->nPixelsCamera; i++){
		branchNameG = TString::Format("Pix[%d]",i);
		treeOutGoldPlated->Branch(branchNameG,
			"GoldPlated",&goldPlatedOut[i],64000,0);
	}

	for(int i = 0; i<readConf->nPixelsCamera; i++){
		branchNameET = TString::Format("Pix[%d]",i);
		treeOutTest->Branch(branchNameET,
			"ExtractedData",&extractedDataTestOut[i],64000,0);
	}
	for(int i = 0; i<readConf->nPixelsCamera; i++){
		branchNameGT = TString::Format("Pix[%d]",i);
		treeOutGoldPlatedTest->Branch(branchNameGT,
			"GoldPlated",&goldPlatedOutTest[i],64000,0);
	}

	

	//cout << "Here 4" <<endl;

	treeOutBiFocal->Branch("TriggerRegion", &bifocalInfoBifocal);
	treeOutGoldPlated->Branch("TriggerRegion", &bifocalInfoGoldPlated);

	treeOutTest->Branch("TriggerRegion",&bifocalInfoTest);
	treeOutGoldPlatedTest->Branch("TriggerRegion",&bifocalInfoGoldPlatedTest);

	treeOutBiFocal->Branch("TriggerTime", &triggerTimeBifocalOut);
	treeOutGoldPlated->Branch("TriggerTime", &triggerTimeGoldPlatedOut);

	treeOutTest->Branch("TriggerTime",&triggerTimeTestOut);
	treeOutTest->Branch("TriggerTime",&triggerTimeGoldPlatedTestOut);
	treeOutCorrectionFactor->Branch("CorrectionFactors",&correctionFactorSave);
	//cout<<"There"<<endl;

	LoadCorrectionFactors(readConf, dataFilePrefix, correctionFactors);
	//fileDownloadReady->cd();

	correctionFactorSave = &correctionFactors;


	//fileDownloadReady->cd();

	
	float pxCorrFactor;

	int maxPxIDPair;
	
	int maxChg = 0;
	int maxPxID = 0;

	for(int i = 0; i<nEventsBiFocal; i++){
		//cout<<"Processing Event: "<<i<<endl;
		
		//fileExtractedData->cd();
		treeInBiFocal->GetEntry(i);

		pixIDsROI = bifocalInfo->GetROIPixelIDs();
		musicIDsROI = bifocalInfo->GetTrigMUSICIDs();

		for(int j = musicIDsROI[0]*8; j<musicIDsROI[0]*8+8; j++){
			if(extractedDataBifocal[j]->GetCharge()/correctionFactors[j]>maxChg){
				maxPxID = j;
				maxChg = extractedDataBifocal[j]->GetCharge()/(correctionFactors[j]);
			}
		}
		//cout<<maxChg<<endl;
		maxChg = 0;
		for(int j = musicIDsROI[0]*8+8; j<musicIDsROI[0]*8+8+8; j++){
			if(extractedDataBifocal[j]->GetCharge()/(correctionFactors[j])>maxChg){
				maxPxIDPair = j;
				maxChg = extractedDataBifocal[j]->GetCharge()/(correctionFactors[j]);
			}
		}

		//cout<<maxChg<<endl;
		if(CheckTriggerGeometry(extractedDataBifocal[maxPxID]->GetTimePeak(),extractedDataBifocal[maxPxID+8]->GetTimePeak())){
			if(CheckIsAboveThreshold(extractedDataBifocal[maxPxID]->GetAmplitude()/(correctionFactors[maxPxID]), 
									extractedDataBifocal[maxPxID+8]->GetAmplitude()/(correctionFactors[maxPxID+8]), 
									readConf->softwareThresholdLvl)){
					if(CheckIsAboveThreshold(extractedDataBifocal[maxPxID]->GetAmplitude()/(correctionFactors[maxPxID]),
											extractedDataBifocal[maxPxIDPair]->GetAmplitude()/(correctionFactors[maxPxIDPair]), 
											readConf->goldPlatedThresholdLvl)) {
						//cout<<"Loading Traces"<<endl;
						LoadTracesForGoldPlated(readConf, i, goldPlatedOut,dataFilePrefix);
							// if(j==0){
							// 	cout<<"Loaded one Trace"<<endl;
							// }
							//cout<<"Event: "<<i<<"is Gold Plated!"<<endl;
						for(int j = 0; j<readConf->nPixelsCamera; j++){
							CopyDataFromExtractedBiFocal(goldPlatedOut[j], extractedDataBifocal[j]);
						}
						//cout<<"Done Loading Traces"<<endl;

						bifocalInfoGoldPlated->Copy(bifocalInfo);
						//bifocalInfoGoldPlated->SetTrigMUSICIDs(bifocalInfo->GetTrigMUSICIDs());
						triggerTimeGoldPlatedOut = triggerTimeBifocal;

						fileDownloadReady->cd();
						treeOutGoldPlated->Fill();
						//cout<<"Event # "<< i <<" in BiFocal Tree is GoldPlated"<<endl; 
			
					}else{
						//cout<<"Event: "<<i<<" is Vanilla"<<endl;
						for(int j = 0; j<readConf->nPixelsCamera; j++){
							CopyDataFromExtractedBiFocal(extractedDataBifocalOut[j], extractedDataBifocal[j]);
						}
						//cout<<extractedDataBifocalOut[0]->GetPedestal()<<endl;
						//cout<<"Done Copying"<<endl;
						fileDownloadReady->cd();
						bifocalInfoBifocal->Copy(bifocalInfo);
						//bifocalInfoBifocal>SetTrigMUSICIDs(bifocalInfo->GetTrigMUSICIDs());
						triggerTimeBifocalOut = triggerTimeBifocal;

						fileDownloadReady->cd();
						treeOutBiFocal->Fill();
					}
			}
			else{
				//cout<<"Event # "<< i <<" in BiFocal Tree is lower than Software Threshold"<<endl; 
			}
		} else{
			//cout<<"Event # "<< i <<" in BiFocal Tree is an accidental"<<endl;
		}

	}

	//cout<<"Next"<<endl;
	for(int i = 0; i<nEventsTest;i++){
		fileExtractedData->cd();
		treeInTest->GetEntry(i);

		pixIDsROI = bifocalInfoTest->GetROIPixelIDs();
		musicIDsROI = bifocalInfoTest->GetTrigMUSICIDs();

		for(int j = musicIDsROI[0]*8; j<musicIDsROI[0]*8+8; j++){
			if(extractedDataTest[j]->GetCharge()*(correctionFactors[j])>maxChg){
				maxPxID = j;
				maxChg = extractedDataTest[j]->GetCharge()*(correctionFactors[j]);
			}
		}
		maxChg = 0;
		for(int j = musicIDsROI[0]*8+8; j<musicIDsROI[0]*8+8+8; j++){
			if(extractedDataBifocal[j]->GetCharge()*(correctionFactors[j])>maxChg){
				maxPxIDPair = j;
				maxChg = extractedDataBifocal[j]->GetCharge()*(correctionFactors[j]);
			}
		}

		if(CheckIsAboveThreshold(extractedDataTest[maxPxID]->GetAmplitude()*(correctionFactors[maxPxID]/correctionFactors[maxPxID]), 
									extractedDataTest[maxPxID+8]->GetAmplitude()*(correctionFactors[maxPxID+8]/correctionFactors[maxPxID+8]), 
									readConf->softwareThresholdLvl)){
			if(CheckIsAboveThreshold(extractedDataTest[maxPxID]->GetAmplitude()*(correctionFactors[maxPxID]/correctionFactors[maxPxID]),
										extractedDataTest[maxPxID+8]->GetAmplitude()*(correctionFactors[maxPxIDPair]/correctionFactors[maxPxIDPair]), 
										readConf->goldPlatedThresholdLvl)) {
						//cout<<"Loading Traces"<<endl;
						LoadTracesForGoldPlated(readConf, i, goldPlatedOutTest,dataFilePrefix);
							// if(j==0){
							// 	cout<<"Loaded one Trace"<<endl;
							// }
						for(int j = 0; j<readConf->nPixelsCamera; j++){
							CopyDataFromExtractedBiFocal(goldPlatedOutTest[j], extractedDataTest[j]);
						}
						bifocalInfoGoldPlatedTest->Copy(bifocalInfoTest);
						triggerTimeGoldPlatedTestOut = triggerTimeTest;

						fileDownloadReady->cd();
						treeOutGoldPlatedTest->Fill();

					}else{
						for(int j = 0; j<readConf->nPixelsCamera; j++){
							CopyDataFromExtractedBiFocal(extractedDataTestOut[j], extractedDataTest[j]);
						}
						//cout<<extractedDataBifocalOut[0]->GetPedestal()<<endl;
						//cout<<"Done Copying"<<endl;
						fileDownloadReady->cd();
						bifocalInfoTestOut->Copy(bifocalInfoTest);
						//bifocalInfoBifocal>SetTrigMUSICIDs(bifocalInfo->GetTrigMUSICIDs());
						triggerTimeTestOut = triggerTimeTest;

						fileDownloadReady->cd();
						treeOutTest->Fill();

					}
				}
			}

	



	fileDownloadReady->cd();
	treeOutHLED = treeInHLED->CloneTree();
	//cout<<"Cloned"<<endl;
	treeOutForced = treeInForced->CloneTree();
	//cout<<"Cloned"<<endl;
	treeOutHLED->Fill();
	treeOutForced->Fill();
	treeOutCorrectionFactor->Fill();
	fileDownloadReady->Write();
	//cout<<"Write"<<endl;
	fileDownloadReady->Close();
	fileExtractedData->Close();
	//f_Correction->Close();


}

void SaveCalibrationCoefficients(ReadConfiguration *readConf, vector<float> correctionFactors, string dataFilePrefix){
	string corrFactorFilename = dataFilePrefix + "_CorrectionFactor.txt";

	ofstream correctionFactorFile;
	correctionFactorFile.open(corrFactorFilename.c_str());

	for(int i = 0; i<readConf->nPixelsCamera; i++){
		correctionFactorFile << correctionFactors[i]<<endl;
		//treeCorrectionFactor->Branch(TString::Format("Pixel[%d]",i),&corrFactor[i]);
	}
	
	correctionFactorFile.close();
}

string SearchForLatestCalibrationFile(string dataFilePrefix){
	string path = dataFilePrefix.substr(0,dataFilePrefix.find_last_of("/\\"));
	struct stat t_stat;
	struct dirent *fileInDir;
	struct tm *timeInfo;
	DIR *directory;
	string fileName, fullFilePath, mostRecentFile;

	mostRecentFile = "";

	int daysPast, hoursPast, minutesPast;
	daysPast = 0;
	hoursPast = 0;
	minutesPast = 0;

	directory = opendir(path.c_str());
	if(directory){
		while((fileInDir = readdir(directory)) != NULL){
			fileName = fileInDir->d_name;
			if(fileName.find("CorrectionFactor")<fileName.length()){
				fullFilePath = path+"/"+fileName;
				stat(fullFilePath.c_str(), &t_stat);
				timeInfo = localtime(&t_stat.st_ctime);

				if(timeInfo->tm_yday >= daysPast && timeInfo->tm_hour >= hoursPast && timeInfo->tm_min >= minutesPast){
					daysPast = timeInfo->tm_yday;
					hoursPast = timeInfo->tm_hour;
					minutesPast = timeInfo->tm_min;

					mostRecentFile = fullFilePath;
				}
			}
		}
	}
	cout<<mostRecentFile<<endl;
	return mostRecentFile;
}

vector<float> CompareCalibrationCoefficients(ReadConfiguration *readConf, vector<float> correctionFactors, string calibrationFile){
	//string correctionFactorFileName = dataFilePrefix + "_CorrectionFactor.txt";
	
	string line;
	ifstream fileCorrectionFactor(calibrationFile.c_str());
	vector<float> correctionFactorsPrev = vector<float>(readConf->nPixelsCamera);
	int i = 0;
	while(getline(fileCorrectionFactor,line)){
		correctionFactorsPrev[i] = stof(line);
		i++;
	}

	fileCorrectionFactor.close();


	float tolerance = 0.05;

	for(int i = 0; i<readConf->nPixelsCamera; i++){
		if(TMath::Abs(correctionFactors[i] - correctionFactorsPrev[i])/correctionFactorsPrev[i] > tolerance){
			correctionFactors[i] = correctionFactorsPrev[i];
		}
	}
	
	return correctionFactors;
}

void HLEDCalibration(ReadConfiguration *readConf, string dataFilePrefix){
	/* Software calibration of SiPMs.
	The response of SiPMs biased to similar bias Voltage should be the same
	with the only variation that arising from Poisson Statistics. 
	The HLED has 5 different amplitudes. 
	*/
	
	string dataFileName = dataFilePrefix + "_Extracted.root";
	//cout<<dataFileName<<endl;
	string mostRecentCalibrationFile;
	TFile *fileIn = new TFile(dataFileName.c_str(),"READ");
	TTree *treeInHLED = (TTree*)fileIn->Get("HLED");
	TTree *treeInForced = (TTree*)fileIn->Get("Forced");
	TTree *treeInBiFocal = (TTree*)fileIn->Get("BiFocal");

	Int_t nEventsHLED, nEventsForced, nEventsBiFocal;
	nEventsHLED = 0;
	nEventsForced = 0;
	nEventsBiFocal = 0;

	nEventsHLED = treeInHLED->GetEntries();
	nEventsForced = treeInForced->GetEntries();
	nEventsBiFocal = treeInBiFocal->GetEntries();

	ExtractedData *extractedDataHLED[readConf->nPixelsCamera];

	vector<float> corrFactor = vector<float>(readConf->nPixelsCamera, 1.0);

	//cout<<"Where"<<endl;

	for (int i = 0; i<readConf->nPixelsCamera; i++){
		corrFactor[i] = 0.0;
		extractedDataHLED[i] = 0;
		extractedDataHLED[i] = new ExtractedData();
		treeInHLED->SetBranchAddress(TString::Format("Pixel[%d]",i),&extractedDataHLED[i]);
	}
	

	float maxChg = 0;


	cout<<"HLED Events"<<nEventsHLED<<endl;
	if(nEventsHLED>0){
		for (int i = 0; i<nEventsHLED; i++){
			treeInHLED->GetEntry(i);
			for(int j = 0; j<readConf->nPixelsCamera; j++){
				corrFactor[j] += extractedDataHLED[j]->GetAmplitude();
			}
		}
		for(int i = 0; i<readConf->nPixelsCamera; i++){
			if(corrFactor[i]>maxChg){
				maxChg = corrFactor[i];
			}
		}
		for(int i = 0; i<readConf->nPixelsCamera; i++){
			corrFactor[i] = corrFactor[i]/maxChg;
			/*if(corrFactor[i]>1.25||corrFactor[i]<0.75){
				corrFactor[i] = 1;
			}*/
			//cout<<"Pixel["<<i<<"] "<<corrFactor[i]<<endl;
		}
	}else{
		for(int i = 0; i<readConf->nPixelsCamera; i++){
			corrFactor[i] = 1.0;
		}
	}
	

	string corrFactorOutName = dataFilePrefix +"_CorrectionFactor.txt";
	if(nEventsHLED<2000 && nEventsHLED != 0){
		cout<<"Looking in"<<endl;
		mostRecentCalibrationFile = SearchForLatestCalibrationFile(dataFilePrefix);
		if(mostRecentCalibrationFile != "" && strcmp(mostRecentCalibrationFile.c_str(),corrFactorOutName.c_str() )!= 0){
			cout<<"Here"<<endl;
			corrFactor = CompareCalibrationCoefficients(readConf, corrFactor, mostRecentCalibrationFile);
			//cout<<"There"<<endl;
		}
	}
	cout<<"Saving"<<endl;


	SaveCalibrationCoefficients(readConf, corrFactor, dataFilePrefix);
	

	cout<<"Closing CorrectionFactor File"<<endl;
	fileIn->Close();
	
	
}

// void HLEDCalibration(ReadConfiguration *readConf, string correctionFilePrefix, Display *display){
// 	/* Software calibration of SiPMs.
// 	The response of SiPMs biased to similar bias Voltage should be the same
// 	with the only variation that arising from Poisson Statistics. 
// 	The HLED has 5 different amplitudes. Will choose the 3rd amplitude as it is
// 	the one more isolated from the rest. This will need automation in the future.
// 	This amplitude occurs between 350 and 850 ADC counts FW
// 	*/

// 	string hledFileName;
// 	string corrFactorFilename = correctionFilePrefix + "_CorrectionFactor.root";
// 	int fileCount = 0;
// 	bool isFile = true;
// 	TTree *treeIn;

// 	double hled_lvl [] = {10,30,100,300,500};
// 	double mirrorTransmission = 0.85;
// 	double pdeTemp = 0.5;
// 	double phAverage, phRMS;
// 	double amplitudeTemp;
	

// 	int nEvents;
// 	//Double_t meanCamera, meanPixel, rmsCamera;
// 	Double_t corrFactor[readConf->nPixelsCamera];

// 	UInt_t signal[readConf->nPixelsCamera][readConf->nSamplesADC];
	
// 	Pulse *pulse;

// 	TProfile *p_pixCorrection = new TProfile("p_pixCorrection","",512,-0.5,511.5,"");
// 	//TH1 *p_val[5];
// 	TProfile *p_val = new TProfile("p_val", "", 512, -0.5, 511.5,"");
// 	double maxValue = 0;
// 	double eventAmplitudes[readConf->nPixelsCamera];
// 	TH1 *h_pixel[readConf->nPixelsCamera];


// 	for (int i = 0; i<readConf->nPixelsCamera; i++){
// 		corrFactor[i] = 1.0;
// 		h_pixel[i] = new TH1F(TString::Format("p_val[%d]",i),"",55,-0.05,1.05);
// 	}
// 	/*for(int i = 0; i < 5; i++){
// 		//p_val[i] = new TH1F(TString::Format("p_val[%d]",i),"",512,-0.5,511.5);
// 		//maxValue[i] = 0;
// 	}*/
	
// 	hledFileName = correctionFilePrefix + TString::Format("_%.4d.root",fileCount);
// 	Double_t norm;
// 	while (isFile){
// 		cout<<"Getting File"<<endl;
		
// 		try{
// 			TFile *fileIn = new TFile(hledFileName.c_str(),"READ");
// 			treeIn = (TTree*)fileIn->Get("Hled_tree");
// 			treeIn->SetBranchAddress("SignalValue",signal);

// 			nEvents = treeIn->GetEntries();

// 			for (int i = 0; i < nEvents; i++){
// 				treeIn->GetEntry(i);

// 				if (i%100 == 0){
// 					cout<<"Getting Event: "<<i<<endl;
// 				}

// 				for (int j = 0; j < readConf->nPixelsCamera; j++){
// 					pulse = new Pulse(readConf, signal[j]);
// 					if (j==0){
// 						norm = pulse->GetCharge();
// 					}
// 					eventAmplitudes[j]=pulse->GetAmplitude();
// 					if(maxValue < pulse->GetAmplitude()){
// 						maxValue = pulse->GetAmplitude();
// 					}
// 					/*h_camDist->Fill(pulse->GetAmplitude(),1.0/512);
// 					h_pixDist[j]->Fill(pulse->GetAmplitude());*/
// 					//p_pixCorrection->Fill(j,pulse->GetCharge()/norm);
// 					//p_val[i%5]->Fill(j,pulse->GetAmplitude());
// 					//h_pixel[j]->Fill(pulse->GetAmplitude());
// 					//corrFactor[j] += (Double_t)pulse->GetCharge()/(norm*(Double_t)nEvents);
// 				}

// 				//maxValue = 1.0;

// 				if(maxValue>300){
// 					for (int j = 0; j < readConf->nPixelsCamera; j++){
// 						eventAmplitudes[j] = eventAmplitudes[j]/maxValue;
// 						p_val->Fill(j,eventAmplitudes[j],1.0);
// 						h_pixel[j]->Fill(eventAmplitudes[j]);
// 						corrFactor[j] = eventAmplitudes[j];
// 					}
// 				}

				

// 				maxValue = 0;
// 			}

// 			fileCount++;

// 			throw exception();
// 		}
// 		catch(exception &e){
// 			cout<<"No more files found"<<endl;
// 			isFile = false;
// 		}	
// 	}

	
// 	TFile *fileOut = new TFile	(corrFactorFilename.c_str(),"RECREATE");
// 	TTree *treeOut = new TTree("t1","Correction Factors");
// 	treeOut->Branch(TString::Format("Correction"),corrFactor,"Correction[512]/D");
	
// 	if(readConf->isDisplay){
// 		display->ShowPlotHLEDCoefficients(p_val);
// 		display->ShowPlotHLEDCoefficientsDistribution(h_pixel,75);
// 	}

	

// 	fileOut->cd();
// 	treeOut->Fill();
// 	//treeOut->Print();
// 	fileOut->Write();
// 	fileOut->Close();

// 	cout<<"Done"<<endl;
	
	
// 	//Basis for a debug and/or Demo flag

// 	/*TCanvas *c = new TCanvas("c","Canvas",3200,2000);

// 	p_pixCorrection->Draw("P");

// 	c->Modified();
// 	c->Update();

// 	c->SaveAs("Histograms.jpg");*/

	
// }

bool CheckFileExistsTString(TString fileName){
	bool isFileExists;
	try{
		TFile *f = new TFile(fileName,"READ");
		throw exception();
		f->Close();
	}
	catch(exception &e){
		cout<<"No such file: "<<fileName<<endl;
		isFileExists = false;
	}

	return isFileExists;
}

void EventPrioritizer(std::string dataFilePrefix, string outputFileDirectory){
	string downloadFileName = dataFilePrefix + "_Download"+".root";

	/*
	used for trimming the date to the right format
	TDatime returns YYYYMMDD which for this year is
	2022MMDD and my getting the modulus of this number
	and 20000000 obtain 22MMDD. Which is the desired format
	*/
	int moduloTrimmer = 20000000;
	
	TDatime *timeOfFile = new TDatime();
	Int_t dateInt, timeInt, priority, index;

	dateInt = timeOfFile->GetDate()%moduloTrimmer;
	timeInt = timeOfFile->GetTime();
	priority = 0;
	index = 0;
	bool isFileExists;

	//cout<<outputFileDirectory<<endl;

	TFile *f_data = new TFile(downloadFileName.c_str(),"READ");
	TTree *treeInBiFocal, *treeInHLED, *treeInForced, *treeInGoldPlated, *treeInCorrectionFactor;
	TTree *treeInTest, *treeInGoldPlatedTest;

	Int_t nEntriesGoldPlated, nEntriesBiFocal, nEntriesHLED, nEntriesTestGoldPlated, nEntriesTest;

	treeInBiFocal = (TTree*)f_data->Get("BiFocal");
	treeInGoldPlated = (TTree*)f_data->Get("GoldPlated");
	treeInTest = (TTree*)f_data->Get("Test");
	treeInGoldPlatedTest = (TTree*)f_data->Get("GoldPlatedTest");
	treeInHLED = (TTree*)f_data->Get("HLED");
	treeInForced= (TTree*)f_data->Get("Forced");
	treeInCorrectionFactor = (TTree*)f_data->Get("CorrectionFactor");

	nEntriesGoldPlated = treeInGoldPlated->GetEntries();
	nEntriesTestGoldPlated = treeInGoldPlatedTest->GetEntries();
	nEntriesTest = treeInTest->GetEntries();
	nEntriesBiFocal = treeInBiFocal->GetEntries();
	nEntriesHLED = treeInHLED->GetEntries();

	if(nEntriesGoldPlated > 0 || nEntriesTestGoldPlated > 0){
		priority = 1;
	}

	TString outFileName;

	do{
		outFileName = TString::Format("%sCT_%d_%d_%.6d_%.5d.dat",outputFileDirectory.c_str(),priority,dateInt,timeInt,index);
		isFileExists = CheckFileExists((string)outFileName.Data());
		index++;
	}while(isFileExists);
	index--;

	//cout<<"Found file"<<endl;
	if(nEntriesGoldPlated > 0){
		cout<<"Saved Download Ready File: "<<outFileName<<endl;
	
		TFile *fileGoldPlatedOut = new TFile(outFileName,"RECREATE");
		TTree *treeOutGoldPlated = treeInGoldPlated->CloneTree();
		TTree *treeOutCorrectionFactorGP = treeInCorrectionFactor->CloneTree();
		treeOutGoldPlated->Fill();
		treeOutCorrectionFactorGP->Fill();
	
		fileGoldPlatedOut->Write();
		fileGoldPlatedOut->Close();

	}

	do{
		outFileName = TString::Format("%sCT_%d_%d_%.6d_%.5d.dat",outputFileDirectory.c_str(),priority,dateInt,timeInt,index);
		isFileExists = CheckFileExists((string)outFileName.Data());
		index++;
	}while(isFileExists);
	index--;

	if(nEntriesTestGoldPlated > 0){
		cout<<"Saved Download Ready File: "<<outFileName<<endl;
	
		TFile *fileGoldPlatedTestOut = new TFile(outFileName,"RECREATE");
		TTree *treeOutGoldPlatedTest = treeInGoldPlatedTest->CloneTree();
		TTree *treeOutCorrectionFactorGPTest = treeInCorrectionFactor->CloneTree();
		treeOutGoldPlatedTest->Fill();
		treeOutCorrectionFactorGPTest->Fill();
	
		fileGoldPlatedTestOut->Write();
		fileGoldPlatedTestOut->Close();

	}

	//cout<<"GoldPlates"<<endl;
	TString outVanillaFileName;
	if(nEntriesHLED > 1000){
		priority = 1;
	}else if(nEntriesBiFocal < 1200){
		priority = 2;
	}else{
		priority = 3;
	}

	//cout<<"Priority: "<<priority<<endl;

	index = 0;
	isFileExists = false;
	do{
		outVanillaFileName = TString::Format("%sCT_%d_%d_%.6d_%.5d.dat",outputFileDirectory.c_str(),priority,dateInt,timeInt,index);
		isFileExists = CheckFileExists((string)outVanillaFileName.Data());
		index++;
	}while(isFileExists);
	index --;

	TFile *fileOutDownload = new TFile(outVanillaFileName,"RECREATE");

	TTree *treeOutHLED, *treeOutBiFocal, *treeOutForced, *treeOutCorrectionFactor;
	//cout<<"Ready for last"<<endl;
	if(nEntriesHLED >= 1000){
		treeOutHLED = treeInHLED->CloneTree();
	}else{
		treeOutHLED = treeInHLED->CloneTree(5);
	}

	cout<<"Saved Download Ready File: "<< outVanillaFileName <<endl;

	treeOutBiFocal = treeInBiFocal->CloneTree();
	treeOutForced = treeInForced->CloneTree();
	treeOutCorrectionFactor = treeInCorrectionFactor->CloneTree();

	fileOutDownload->cd();

	treeOutHLED->Fill();
	treeOutBiFocal->Fill();
	treeOutForced->Fill();
	treeOutCorrectionFactor->Fill();
	
	fileOutDownload->Write();
	fileOutDownload->Close();

	
	f_data->cd();
	f_data->Close();

	//cout<<"Bye"<<endl;
}


/*void CombineFiles(ReadConfiguration *readConf, string dataFilePrefix, string correctionFilePrefix){
	string extractedFileName = dataFilePrefix + "_Extracted"+".root";
	string corrFactorFilename = correctionFilePrefix + "_CorrectionFactor.root";


	
	used for trimming the date to the right format
	TDatime returns YYYYMMDD which for this year is
	2022MMDD and my getting the modulus of this number
	and 20000000 obtain 22MMDD. Which is the desired format
	
	int moduloTrimmer = 20000000;
	
	TDatime *timeOfFile = new TDatime();
	Int_t dateInt, timeInt, priority, index;

	dateInt = timeOfFile->GetDate()%moduloTrimmer;
	timeInt = timeOfFile->GetTime();
	priority = 0;
	index = 0;

	bool isFileExists;

	TString outFileName;

	do{
		outFileName = TString::Format("CT_%d_%d_%d_%.5d.bin",priority,dateInt,timeInt,index);
		isFileExists = CheckFileExists((string)outFileName.Data());
		index++;
	}while(isFileExists);

	ExtractedData *ex[readConf->nPixelsCamera], *ex_out[readConf->nPixelsCamera];
	Double_t corrFactor[readConf->nPixelsCamera];

	TFile *f_data = new TFile(extractedFileName.c_str(),"READ");
	TTree *t_data = (TTree*)f_data->Get("t1");


	TFile *f_corr = new TFile(corrFactorFilename.c_str(),"READ");
	TTree *t_corr = (TTree*)f_corr->Get("t1");


	TFile *f_out = new TFile(outFileName,"RECREATE");
	TTree *t_out = new TTree("t1","Merged Data");



	for(int i = 0; i<readConf->nPixelsCamera; i++){
		ex_out[i] = 0;
		ex[i] = 0;
		t_out->Branch(TString::Format("Pixel[%d]",i),
			"ExtractedData",&ex_out[i],64000,0);
		t_data->SetBranchAddress(TString::Format("Pixel[%d]",i),
			&ex[i]);
	}
	
	t_corr->SetBranchAddress("Correction",&corrFactor);
	t_corr->GetEntry(0);

	for(int i = 0; i<t_data->GetEntries(); i++){
		t_data->GetEntry(i);
		for(int j = 0; j<readConf->nPixelsCamera; j++){
			ex_out[j]->SetAmplitude(ex[j]->GetAmplitude());
			ex_out[j]->SetTimePeak(ex[j]->GetTimePeak());
			ex_out[j]->SetTimeFWHM(ex[j]->GetTimeFWHM());
			ex_out[j]->SetTimeFW(ex[j]->GetTimeFW());
			ex_out[j]->SetCorrFactor(corrFactor[j]);

		}
		f_out->cd();
		t_out->Fill();
	}
	f_out->cd();
	f_out->Write();
	f_out->Close();


}*/

// void OutFilesDownload(ReadConfiguration *readConf, string dataFilePrefix, string outputFileDirectory){
// 	string extractedFileName = dataFilePrefix + "_Extracted"+".root";
// 	//string corrFactorFilename = correctionFilePrefix + "_CorrectionFactor.root";


// 	/*
// 	used for trimming the date to the right format
// 	TDatime returns YYYYMMDD which for this year is
// 	2022MMDD and my getting the modulus of this number
// 	and 20000000 obtain 22MMDD. Which is the desired format
// 	*/
// 	int moduloTrimmer = 20000000;
	
// 	TDatime *timeOfFile = new TDatime();
// 	Int_t dateInt, timeInt, priority, index;

// 	dateInt = timeOfFile->GetDate()%moduloTrimmer;
// 	timeInt = timeOfFile->GetTime();
// 	priority = 0;
// 	index = 0;
// 	bool isFileExists;

// 	//cout<<outputFileDirectory<<endl;
// 	TString outFileName;

// 	do{
// 		outFileName = TString::Format("%sCT_%d_%d_%.6d_%.5d.dat",outputFileDirectory.c_str(),priority,dateInt,timeInt,index);
// 		isFileExists = CheckFileExists((string)outFileName.Data());
// 		index++;
// 	}while(isFileExists);

// 	ExtractedData *ex[readConf->nPixelsCamera], *ex_out[readConf->nPixelsCamera];
	
// 	TFile *f_data = new TFile(extractedFileName.c_str(),"READ");
// 	TTree *t_data = (TTree*)f_data->Get("t1");


// 	TFile *f_out = new TFile(outFileName,"RECREATE");
// 	TTree *t_out = new TTree("t1","Merged Data");



// 	for(int i = 0; i<readConf->nPixelsCamera; i++){
// 		ex_out[i] = 0;
// 		ex[i] = 0;
// 		t_out->Branch(TString::Format("Pixel[%d]",i),
// 			"ExtractedData",&ex_out[i],64000,0);
// 		t_data->SetBranchAddress(TString::Format("Pixel[%d]",i), &ex[i]);
// 	}
	

// 	for(int i = 0; i<t_data->GetEntries(); i++){
// 		t_data->GetEntry(i);
// 		for(int j = 0; j<readConf->nPixelsCamera; j++){
// 			ex_out[j]->SetAmplitude(ex[j]->GetAmplitude());
// 			ex_out[j]->SetTimePeak(ex[j]->GetTimePeak());
// 			/*ex_out[j]->SetTimeFWHM(ex[j]->GetTimeFWHM());
// 			ex_out[j]->SetTimeFW(ex[j]->GetTimeFW());*/
// 		}

// 		f_out->cd();
// 		t_out->Fill();
// 	}
// 	f_out->cd();
// 	f_out->Write();
// 	f_out->Close();
// }

void PrintHelpMessage(){
	 cout << "ExACT" << endl;
     cout << endl;
     cout << "Command Line options: " << endl;
     cout << "\t -h or --help:\t\t\tThis help message"<<endl;
     cout << "\t -df or --datafile <Data File>\t\t\tInput file coming from EventBuilder"<<endl;
     cout << "\t -c or --configfile <Config File>\t\t\tConfiguration file ExACT.cfg"<<endl;
     cout << "\t -od or --outputdir\t\t\tFolder where to store the files generated by ExACT"<<endl;
     cout << endl;
     exit( 0 );
}

bool StringToBool(string sBool){
	bool isValueBool;
	transform(sBool.begin(), sBool.end(), sBool.begin(), ::toupper);
	if(sBool=="1" || sBool=="TRUE"){
		isValueBool = true;
	}
	else{
		isValueBool = false;
	}

	return isValueBool;
}
int main (int argc, char **argv){
	auto start = chrono::high_resolution_clock::now();

	string flag;
	string dataFilePrefix, cfgFileName, outputFileDirectory;
	bool isUnderWater;
	ReadConfiguration *readConf;
	string underWater;

	if(argc<3){
		PrintHelpMessage();
	}else{
		for(int i = 1; i<argc; i++){
			flag = argv[i];
			if(i+1>argc){
				cout<<"Argument for flag " << argv[i] << " is not present"<<endl;
			}else{
				if(flag=="-h" || flag == "--help"){
					PrintHelpMessage();
				}else if(flag=="-df" || flag == "--datafile"){
					dataFilePrefix = argv[i+1];
				}else if(flag=="-c" || flag == "--configfile"){
					cfgFileName = argv[i+1];
					readConf = new ReadConfiguration(cfgFileName);
				}else if(flag=="-od" || flag == "--outputdir"){
					outputFileDirectory = argv[i+1];
				}else if(flag=="-uw" || flag == "--offline"){
					underWater = argv[i+1];
					isUnderWater = StringToBool(underWater);
				}
			}
			i++;
		}
	}
	
	cout << dataFilePrefix << " " << cfgFileName << " " <<outputFileDirectory<<endl;
	
	if(!isUnderWater){
		
		// Adding a Starlink file preparation. It will keep traces for all events
		// that fulfill bifocal condition.
		StarlinkFileExtraction(readConf, dataFilePrefix);

		// First Extract all fo the events obtained for the data arquisition
		// run. Every event is reduced to charge and amplitude which will be
		// calibrated, discriminated and prioritized in later steps.

		PulseExtraction(readConf,dataFilePrefix);

		//Prepare HLED correction Factors for this run.

		/*If HLED calibration finds only HLED events and # events exceeds
		400 per amplitude. Will consider it an HLED calibration run and will 
		generate a file with correction factors*/

		HLEDCalibration(readConf, dataFilePrefix);

		//Discriminator

		// Look at BiFocal events and remove accidental triggers from 
		// dataset. Algorithm to looks at the triggered MUSICs,
		// and find the maximum amplitude signal. Check if there is a corresponding
		// signal in the pixel expected.
		// Impose a software threshold to ensure that both signals are a copy of 
		// each other.

		EventDiscrimination(readConf, dataFilePrefix);

		//Prioritizer

		EventPrioritizer(dataFilePrefix, outputFileDirectory);

		//File for download
		// OutFilesDownload(readConf,dataFilePrefix,outputFileDirectory);
	}else{
		// Here is where you will write your own modular methods //
		/* IPlotTools.h and IPlotTools.cpp show how modular methods for analysis should be 
		implemented.

		We want each stage or analysis part to be run in here as 

		AnalysisTool::ObtainInformation(foo, bar)

		To do so, a .h (header) file needs to be created in the include directory
		The definitions will go in the .cpp file in the src directory

		EXAMPLE:
		
		--------------
		AnalysisTool.h
		--------------

		class AnalysisTool{
			public:
				static int ObtainInformation(int foo, int bar);
		};

		-----------------
		AnalysisTool.cpp
		-----------------

		int ObtainInformation(int foo, int bar){
			return foo*bar;
		}
		

		By doing it this way, the dependency of different modules with each other are decreased.
		*/
	}


	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop-start).count();
	duration *= 1e-3;

	cout<<"Time Elapsed: "<<duration<<" seconds"<<endl;


	exit(0);
}