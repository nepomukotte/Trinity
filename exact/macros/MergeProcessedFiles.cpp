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
          found=file.find(".dat");
          if (found==string::npos) continue;
          names->push_back(file);
        }
        gSystem->FreeDirectory(dir);
    }
     // for(size_t i(0); i<names->size(); i++)
     //   Printf("%s",(*names)[i].c_str());
}

int MergeProcessedFiles (std::string MainDir,std::string fileNameOut){

    string DataDir = MainDir;
    vector<string> * filenames = new vector<string>;
    readDir(DataDir,filenames);
    

    string outFile =  MainDir+fileNameOut;

    TFile *fl_merge = new TFile(outFile.c_str(),"RECREATE");
    TTree *tHLED = new TTree("HLED","HLED Triggers");
    /*TTree *tBiFocal = new TTree("BiFocal","BiFocal Triggers");
    TTree *tForced = new TTree("Forced","Forced Triggers");
    TTree *tTest = new TTree("Test","Test Triggers");
*/
	ExtractedData *evHLEDOut[512];
	/*ExtractedData *evBiFocalOut[512];
	ExtractedData *evForcedOut[512];
	ExtractedData *evTestOut[512];
*/
    unsigned long long eventTimeIn, eventTimeOut;

    for(int i = 0; i<512; i++){
        evHLEDOut[i] = new ExtractedData();
        /*evBiFocalOut[i] = new ExtractedData();
        evForcedOut[i] = new ExtractedData();
        evTestOut[i] = new ExtractedData();*/
    }


    for(int i = 0; i<512; i++){
        tHLED->Branch(TString::Format("Pixel[%d]",i),"ExtractedData", &evHLEDOut[i],64000,0);
       /* tBiFocal->Branch(TString::Format("Pixel[%d]",i),"ExtractedData", &evBiFocalOut[i],64000,0);
        tForced->Branch(TString::Format("Pixel[%d]",i),"ExtractedData", &evForcedOut[i],64000,0);
        tTest->Branch(TString::Format("Pixel[%d]",i),"ExtractedData", &evTestOut[i],64000,0);*/

    }
	tHLED->Branch("TriggerTime", &eventTimeOut);



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
    /*TTree *tBiFocalIn = 0;
    TTree *tForcedIn = 0;
    TTree *tTestIn = 0;*/
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
	ExtractedData *evHLED[512];
	/*ExtractedData *evBiFocal[512];
	ExtractedData *evForced[512];
	ExtractedData *evTest[512];*/

    

    for(int i = 0; i<512; i++){
        evHLED[i] = 0;
        /*evBiFocal[i] = 0;
        evForced[i] =0;
        evTest[i] = 0;*/
    }

	
    while(source){

    	source->cd();
        tHLEDIn = (TTree*)source->Get("HLED");
       /* tForcedIn = (TTree*)source->Get("Forced");
        tTestIn = (TTree*)source->Get("Test");
        tBiFocalIn = (TTree*)source->Get("BiFocal");*/
        

        for(int i = 0; i<512; i++){
            tHLEDIn->SetBranchAddress(TString::Format("Pixel[%d]",i),&evHLED[i]);
            
            /*tBiFocalIn->SetBranchAddress(TString::Format("Pixel[%d]",i),&evBiFocal  [i]);
            tTestIn->SetBranchAddress(TString::Format("Pixel[%d]",i),&evTest[i]);
            tForcedIn->SetBranchAddress(TString::Format("Pixel[%d]",i),&evForced[i  ]);*/
        }
        tHLEDIn->SetBranchAddress(TString::Format("TriggerTime"),&eventTimeIn);
	
	   for(int i = 0; i<tHLEDIn->GetEntries(); i++){
            tHLEDIn->GetEntry(i);
            for(int j = 0; j<512; j++){
                evHLEDOut[i] = evHLED[i];
            }
            eventTimeOut = eventTimeIn;
	   	
            fl_merge->cd();
            tHLED->Fill();

            if(i%100==0){
                cout<<"Processed Event: "<<i<<endl;
            }
	   }

	/*for(int i = 0; i<tBiFocalIn->GetEntries(); i++){
		tBiFocalIn->GetEntry(i);
        for(int j = 0; j<512; j++){
		  evBiFocalOut[i] = evBiFocal[i];
		}
        fl_merge->cd();
		tBiFocal->Fill();
	}

	for(int i = 0; i<tForcedIn->GetEntries(); i++){
		tForcedIn->GetEntry(i);
        for(int j = 0; j<512; j++){
		  evForcedOut[i] = evForced[i];
        }
		fl_merge->cd();
		tForced->Fill();
	}

	for(int i = 0; i<tTestIn->GetEntries(); i++){
		tTestIn->GetEntry(i);
        for(int j = 0; j<512; j++){
            evTestOut[i] = evTest[i];
        }
		fl_merge->cd();
		tTest->Fill();
	}*/


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
    //fl_merge->Write();
    fl_merge->Close();


    return 0;


}
