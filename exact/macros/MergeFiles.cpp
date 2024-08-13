R__LOAD_LIBRARY(libExACT.so)

#include <TTree.h>
#include <TFile.h>
//#include <Event.h>
void readDir(string dirname,vector<string> *names) {

	string file;
	void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname.c_str()));
    const char* strptr;
    if (dir) {
        while ((strptr = gSystem->GetDirEntry(dir))) {
          file = strptr;
          cout<<file<<endl;
          //if (file=="." || file==".." || file=="log") continue;
          //names->push_back(file);
          size_t found;
          found=file.find(".root");
          if (found==string::npos) continue;
          names->push_back(file);
        }
        gSystem->FreeDirectory(dir);
    }
     // for(size_t i(0); i<names->size(); i++)
     //   Printf("%s",(*names)[i].c_str());
}

int MergeFiles (std::string MainDir,std::string fileNameOut){

    string DataDir = MainDir;
    vector<string> * filenames = new vector<string>;
    readDir(DataDir,filenames);
    

    string outFile =  MainDir+fileNameOut;

    TFile *fl_merge = new TFile(outFile.c_str(),"RECREATE");
    TTree *tHLED = new TTree("HLED","HLED Triggers");
    TTree *tBiFocal = new TTree("BiFocal","BiFocal Triggers");
    TTree *tForced = new TTree("Forced","Forced Triggers");
    TTree *tTest = new TTree("Test","Test Triggers");

	Event *evHLEDOut = new Event();
	Event *evBiFocalOut = new Event();
	Event *evForcedOut = new Event();
	Event *evTestOut = new Event();

	tHLED->Branch("Events","Event", &evHLEDOut,64000,0);
	tBiFocal->Branch("Events","Event", &evBiFocalOut,64000,0);
	tForced->Branch("Events","Event", &evForcedOut,64000,0);
	tTest->Branch("Events","Event", &evTestOut,64000,0);



    /*TChain *ch_Glob = 0;
    TChain *ch_T0_Glob = 0;*/
    //TDirectory *dir_events = fl_merge->mkdir("Events","Events");
    //fl_merge->cd("Events");


    /*int iNumPixels = 512;

    TTree *t_merg_tSim = new TTree("tSimulatedEvents", "tSimulatedEvents");
    TTree *t_merg_T0 = new TTree("T0","T0");

    UInt_t uNumTelescopes;
    Float_t energy;
    UInt_t eventNumber ;
    Float_t xcore ;
    Float_t ycore ;
    Float_t azPrim ;
    Float_t znPrim ;
    Bool_t arrayTriggerBit ;
    std::vector< Bool_t > *vTelescopeTriggerBits =0 ;
    Float_t DeltaTL3 ;

    UInt_t uNumTelescopes_t;
    Float_t energy_t;
    UInt_t eventNumber_t ;
    Float_t xcore_t;
    Float_t ycore_t;
    Float_t azPrim_t;
    Float_t znPrim_t;
    Bool_t arrayTriggerBit_t;
    std::vector< Bool_t > *vTelescopeTriggerBits_t =0 ;
    Float_t DeltaTL3_t;
    vector< vector<Int_t> *>   iFADCTraceInPixel;
    vector<Int_t>   *iPEInPixel;

    vector<int> *vTriggerCluster=0;
    vector<Float_t> *fTimeOverThreshold=0;
    vector<Float_t>   *fSumTimeInPixel=0;
    vector<Int_t>   *iQDCInPixel=0;
    Int_t           iNumPhotonsInFocalPlane=0;
    Float_t         fAzTel=0 ;
    Float_t         fZnTel =0;
    vector<Bool_t>  *bInLoGain=0;
    iFADCTraceInPixel.assign(iNumPixels,0);
    cout<<&iFADCTraceInPixel[0]<<endl;

    vector<int> *vTriggerCluster_t=0;
    vector<Float_t> *fTimeOverThreshold_t=0;
    vector<Float_t>   *fSumTimeInPixel_t=0;
    vector<Int_t>   *iQDCInPixel_t=0;
    Int_t           iNumPhotonsInFocalPlane_t=0;
    Float_t         fAzTel_t=0 ;
    Float_t         fZnTel_t=0;
    vector<Bool_t>  *bInLoGain_t=0;
    vector< vector<Int_t> *>   iFADCTraceInPixel_t;
    vector<Int_t>   *iPEInPixel_t=0;
    iFADCTraceInPixel_t.assign(iNumPixels,0);
    cout<<&iFADCTraceInPixel_t[0]<<endl;

    t_merg_tSim->Branch("energy",&energy);
    t_merg_tSim->Branch("ZnPrim",&znPrim);
    t_merg_tSim->Branch("AzPrim",&azPrim);
    t_merg_tSim->Branch("xcore",&xcore);
    t_merg_tSim->Branch("ycore",&ycore);
    t_merg_tSim->Branch("arrayTriggerBit",&arrayTriggerBit);
    t_merg_tSim->Branch("uNumTelescopes",&uNumTelescopes);
    t_merg_tSim->Branch("eventNumber",&eventNumber);
    t_merg_tSim->Branch("DeltaTL3",&DeltaTL3);
    t_merg_tSim->Branch("vTelescopeTriggerBits",&vTelescopeTriggerBits);

    t_merg_T0->Branch("vGroupsInTriggerCluster",&vTriggerCluster);
    t_merg_T0->Branch("vTimeOverThreshold",&fTimeOverThreshold);
    t_merg_T0->Branch("vSumTimeInPixel", &fSumTimeInPixel);
    t_merg_T0->Branch("vPEInPixel", &iPEInPixel);
    t_merg_T0->Branch("vQDCValue", &iQDCInPixel);
    t_merg_T0->Branch("iPhotonsInFocalPlane", &iNumPhotonsInFocalPlane);
    t_merg_T0->Branch("fAzTel", &fAzTel);
    t_merg_T0->Branch("fZnTel", &fZnTel);
    t_merg_T0->Branch("vHiLoGainBit", &bInLoGain);
    TString name;
    for(int g=0;g<iNumPixels;g++){
            name.Form("vFADCTraces%i",g);
            //cout<<name<<endl;
            t_merg_T0->Branch(name,&iFADCTraceInPixel[g]);
    }*/

    TList *lst_Files = new TList();

    for (int i = 0; i<filenames->size(); i++){
            stringstream fnamestream;
            fnamestream<<DataDir<<"/"<<(*filenames)[i].c_str();
            string fstring(fnamestream.str());
            cout<<fstring<<endl;
            lst_Files->Add(TFile::Open(fstring.c_str()));
    }


    //TString path( (char*)strstr( dir_events->GetPath(), ":" ) );
    //path.Remove( 0, 2 );


//    cout<<"Working Path: "<<path<<endl;


    TFile *source = (TFile*)lst_Files->First();
    source->cd();
    TTree *tHLEDIn= 0;
    TTree *tBiFocalIn = 0;
    TTree *tForcedIn = 0;
    TTree *tTestIn = 0;
//    TTree *t_temp_tSim;
//    TTree *t_temp_T0;

    //fstream fileHLED, fileBiFocal, fileForced, fileTest;
    //std::string fnameHLED, fnameBiFocal, fnameForced,fnameTest;
    //fnameHLED= MainDir+"/MergedHLED.dat";
    //fnameBiFocal= MainDir+"/MergedBiFocal.dat";
    //fnameForced= MainDir+"/MergedForced.dat";
    //fnameTest= MainDir+"/MergedTest.dat";

//    fileHLED.open(fnameHLED.c_str(),ios::out|ios::binary);
//    fileBiFocal.open(fnameBiFocal.c_str(),ios::out|ios::binary);
//    fileForced.open(fnameForced.c_str(),ios::out|ios::binary);
//    fileTest.open(fnameTest.c_str(),ios::out|ios::binary);
	Event *evHLED = 0;
	Event *evBiFocal = 0;
	Event *evForced = 0;
	Event *evTest = 0;

	
    while(source){

    	source->cd();
        tHLEDIn = (TTree*)source->Get("HLED");
        tForcedIn = (TTree*)source->Get("Forced");
        tTestIn = (TTree*)source->Get("Test");
        tBiFocalIn = (TTree*)source->Get("BiFocal");
        

	tHLEDIn->SetBranchAddress("Events",&evHLED);
	tBiFocalIn->SetBranchAddress("Events",&evBiFocal);
	tTestIn->SetBranchAddress("Events",&evTest);
	tForcedIn->SetBranchAddress("Events",&evForced);
	
	for(int i = 0; i<tHLEDIn->GetEntries(); i++){
		tHLEDIn->GetEntry(i);
		evHLEDOut = evHLED;
		fl_merge->cd();
		tHLED->Fill();
	}

	for(int i = 0; i<tBiFocalIn->GetEntries(); i++){
		tBiFocalIn->GetEntry(i);
		evBiFocalOut = evBiFocal;
		fl_merge->cd();
		tBiFocal->Fill();
	}

	for(int i = 0; i<tForcedIn->GetEntries(); i++){
		tForcedIn->GetEntry(i);
		evForcedOut = evForced;
		fl_merge->cd();
		tForced->Fill();
	}

	for(int i = 0; i<tTestIn->GetEntries(); i++){
		tTestIn->GetEntry(i);
		evTestOut = evTest;
		fl_merge->cd();
		tTest->Fill();
	}


	//fileHLED.write((char*)tHLEDIn,sizeof(tHLEDIn));
	//fileForced.write((char*)tForcedIn,sizeof(tForcedIn));
	//fileTest.write((char*)tTestIn,sizeof(tTestIn));
	//fileBiFocal.write((char*)tBiFocalIn,sizeof(tBiFocalIn));

//        tHLED = tHLEDIn->CloneTree();
//        tBiFocal = tBiFocalIn->CloneTree();
//        tForced = tForcedIn->CloneTree();
//        tTest = tTestIn->CloneTree();

/*        fl_merge->cd();
        tHLED->Fill();
        tBiFocal->Fill();
        tTest->Fill();
        tForced->Fill();
*/
        source = (TFile*)lst_Files->After(source);
    }
	
//	fileHLED.close();
//	fileForced.close();
//	fileTest.close();
//	fileBiFocal.close();
    fl_merge->cd();
    fl_merge->Write();
    fl_merge->Close();


    return 0;


}
