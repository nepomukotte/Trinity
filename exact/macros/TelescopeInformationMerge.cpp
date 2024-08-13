R__LOAD_LIBRARY(libExACT.so)

#include <TTree.h>
#include <TFile.h>


int TelescopeInformationMerge (std::string filename1, std::string filename2, std::string filenameTilt, std::string filenameTiltRaw ,std::string filename, std::string fileOut){

	std::vector<std::string> headers;
	std::vector<std::vector<std::string>> data;
	std::vector<std::string> headers2;
	std::vector<std::vector<std::string>> data2;
    std::vector<std::string> headersTilt;
    std::vector<std::vector<std::string>> dataTilt;
    std::vector<std::string> headersTiltRaw;
    std::vector<std::vector<std::string>> dataTiltRaw;

	std::ifstream file1(filename1);

    int timestamp2Index, timestampIndex;
    if (!file1.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename1);
    }

    std::string line;
    if (std::getline(file1, line)) {
        std::stringstream ss(line);
        std::string header;
        while (ss >> header) {
            headers.push_back(header);
        }
    }

    if (headers.empty() || headers.size() != 9) {
        throw std::runtime_error("Invalid file format: " + filename1);
    }

    timestampIndex = 0;
    for (int i = 0; i < headers.size(); ++i) {
        if (headers[i] == "timestamp") {
            timestampIndex = i;
            break;
        }
    }

    while (std::getline(file1, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;

        while (ss >> cell) {
            row.push_back(cell);
        }

        if (row.size() != headers.size()) {
            throw std::runtime_error("Invalid file format: " + filename1);
        }

        data.push_back(row);
    }

    file1.close();

    std::ifstream file2(filename2);
    if (!file2.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename2);
    }

    if (std::getline(file2, line)) {
        std::stringstream ss(line);
        std::string header;
        while (ss >> header) {
            headers2.push_back(header);
        }
    }

    if (headers2.empty()) {
        throw std::runtime_error("Invalid file format: " + filename2);
    }

    timestamp2Index = 0;
    for (int i = 0; i < headers2.size(); ++i) {
        if (headers2[i] == "timestamp2") {
            timestamp2Index = i;
            break;
        }
    }

    while (std::getline(file2, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;

        while (ss >> cell) {
            row.push_back(cell);
        }

        if (row.size() != headers2.size()) {
            throw std::runtime_error("Invalid file format: " + filename2);
        }

        data2.push_back(row);
    }

    file2.close();


    std::ifstream file3(filenameTilt);
    if (!file3.is_open()) {
        throw std::runtime_error("Failed to open file: " + filenameTilt);
    }

    if (std::getline(file3, line)) {
        std::stringstream ss(line);
        std::string header;
        while (ss >> header) {
            headersTilt.push_back(header);
        }
    }

    if (headersTilt.empty()) {
        cout<<"Header"<<endl;
        throw std::runtime_error("Invalid file format: " + filenameTilt);
    }

    timestamp2Index = 0;
   /* for (int i = 0; i < headers2.size(); ++i) {
        if (headers2[i] == "timestamp") {
            timestamp2Index = i;
            break;
        }
    }*/

    while (std::getline(file3, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;

        while (ss >> cell) {
            row.push_back(cell);
        }

        if (row.size() != headersTilt.size()) {
            cout<<"Row"<<endl;
            throw std::runtime_error("Invalid file format: " + filenameTilt);
        }

        dataTilt.push_back(row);
    }

    file3.close();

    std::ifstream file4(filenameTiltRaw);
    if (!file4.is_open()) {
        throw std::runtime_error("Failed to open file: " + filenameTiltRaw);
    }

    if (std::getline(file4, line)) {
        std::stringstream ss(line);
        std::string header;
        while (ss >> header) {
            headersTiltRaw.push_back(header);
        }
    }

    if (headersTilt.empty() || headersTilt.size() != 2) {
        cout<<"Header"<<endl;
        throw std::runtime_error("Invalid file format: " + filenameTiltRaw);
    }

    timestamp2Index = 0;
   /* for (int i = 0; i < headers2.size(); ++i) {
        if (headers2[i] == "timestamp") {
            timestamp2Index = i;
            break;
        }
    }*/

    while (std::getline(file4, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;

        while (ss >> cell) {
            row.push_back(cell);
        }

        if (row.size() != headersTiltRaw.size()) {
            cout<<"Row"<<endl;
            throw std::runtime_error("Invalid file format: " + filenameTiltRaw);
        }

        dataTiltRaw.push_back(row);
    }

    file4.close();

    ULong64_t normalizedTime = 1683936000;
    ULong64_t normalizedTimeNs = 1683936000*1e8; //Should be in 10s of ns
    ULong64_t timeAfterLaunch = 0;
    IEvent *iEvHLED = new IEvent();
    IEvent *iEvBiFocal = new IEvent();
    IEvent *iEvTest = new IEvent();
    IEvent *iEvForced = new IEvent();

    Event *evHLED = 0;
    Event *evBiFocal = 0;
    Event *evTest = 0;
    Event *evForced = 0;


    TFile *fIn = new TFile(filename.c_str(), "READ");
    
    TTree *tHLED = (TTree*) fIn->Get("HLED");    
    TTree *tBiFocal = (TTree*) fIn->Get("BiFocal");    
    TTree *tForced = (TTree*) fIn->Get("Forced");    
    TTree *tTest = (TTree*) fIn->Get("Test");


    tHLED->SetBranchAddress("Events",&evHLED);
    tBiFocal->SetBranchAddress("Events",&evBiFocal);
    tForced->SetBranchAddress("Events",&evForced);
    tTest->SetBranchAddress("Events",&evTest);

	TFile *fOut = new TFile(fileOut.c_str(),"RECREATE");

    TTree *tMHLED = new TTree("HLED","HLED Events");
    TTree *tMBiFocal = new TTree("BiFocal","BiFocal Events");
    TTree *tMForced = new TTree("Forced","Forced Events");
    TTree *tMTest = new TTree("Test","Test Events");
	
    tMHLED->Branch("Events","IEvent",&iEvHLED,64000,99);
    tMBiFocal->Branch("Events","IEvent",&iEvBiFocal,64000,99);
    tMForced->Branch("Events","IEvent",&iEvForced,64000,99);
    tMTest->Branch("Events","IEvent",&iEvTest,64000,99);

    for(int i = 0; i < tHLED->GetEntries(); i++){
        if(i%100 == 0){
            cout << "Processing HLED event: "<<i<<endl;
        }
        tHLED->GetEntry(i);
        timeAfterLaunch = evHLED->GetTBTime() - normalizedTimeNs;
        timeAfterLaunch = timeAfterLaunch/1e8;
        // cout<<evHLED->GetTBTime()<<endl;
        // cout<<normalizedTimeNs<<endl;
        // cout<<timeAfterLaunch<<endl;;
        iEvHLED->SetCoBoTime(evHLED->GetCoBoTime());
        iEvHLED->SetUNIXTime(evHLED->GetUNIXTime());
        iEvHLED->SetTBTime(evHLED->GetTBTime());

        iEvHLED->SetEventType(0);
        iEvHLED->SetROIPixelID(evHLED->GetROIPixelID());
        iEvHLED->SetROIMusicID(evHLED->GetROIMusicID());
        iEvHLED->SetSignalValue(evHLED->GetSignalValue());

        iEvHLED->SetParametersFromTimestamp(timeAfterLaunch,data,data2);
        iEvHLED->SetTelescopePointing(timeAfterLaunch,dataTilt);
        iEvHLED->SetTelescopePointingRaw(timeAfterLaunch,dataTiltRaw);


        tMHLED->Fill();
    }
    for(int i = 0; i < tBiFocal->GetEntries(); i++){
        if(i%100 == 0){
            cout << "Processing BiFocal event: "<<i<<endl;
        }
        fIn->cd();
        tBiFocal->GetEntry(i);
        timeAfterLaunch = evBiFocal->GetTBTime()- normalizedTimeNs;
        timeAfterLaunch = timeAfterLaunch/1e8;
        iEvBiFocal->SetCoBoTime(evBiFocal->GetCoBoTime());
        iEvBiFocal->SetUNIXTime(evBiFocal->GetUNIXTime());
        iEvBiFocal->SetTBTime(evBiFocal->GetTBTime());

        iEvBiFocal->SetEventType(0);
        iEvBiFocal->SetROIPixelID(evBiFocal->GetROIPixelID());
        iEvBiFocal->SetROIMusicID(evBiFocal->GetROIMusicID());
        iEvBiFocal->SetSignalValue(evBiFocal->GetSignalValue());

        iEvBiFocal->SetParametersFromTimestamp(timeAfterLaunch,data,data2);
        iEvBiFocal->SetTelescopePointing(timeAfterLaunch,dataTilt);
        iEvBiFocal->SetTelescopePointingRaw(timeAfterLaunch,dataTiltRaw);

        fOut->cd();
        tMBiFocal->Fill();
    }
    for(int i = 0; i < tForced->GetEntries(); i++){
        if(i%100 == 0){
            cout << "Processing Forced event: "<<i<<endl;
        }
        fIn->cd();
        tForced->GetEntry(i);
        timeAfterLaunch = evForced->GetTBTime()- normalizedTimeNs;
        timeAfterLaunch = timeAfterLaunch/1e8;
        iEvForced->SetCoBoTime(evForced->GetCoBoTime());
        iEvForced->SetUNIXTime(evForced->GetUNIXTime());
        iEvForced->SetTBTime(evForced->GetTBTime());

        iEvForced->SetEventType(0);
        iEvForced->SetROIPixelID(evForced->GetROIPixelID());
        iEvForced->SetROIMusicID(evForced->GetROIMusicID());
        iEvForced->SetSignalValue(evForced->GetSignalValue());

        iEvForced->SetParametersFromTimestamp(timeAfterLaunch,data,data2);
        iEvForced->SetTelescopePointing(timeAfterLaunch,dataTilt);
        iEvForced->SetTelescopePointingRaw(timeAfterLaunch,dataTiltRaw);


        fOut->cd();
        tMForced->Fill();
    }
    for(int i = 0; i < tTest->GetEntries(); i++){
        if(i%100 == 0){
            cout << "Processing Test event: "<<i<<endl;
        }
        fIn->cd();
        tTest->GetEntry(i);
        timeAfterLaunch = evTest->GetTBTime()- normalizedTimeNs;
        timeAfterLaunch = timeAfterLaunch/1e8;
        iEvTest->SetCoBoTime(evTest->GetCoBoTime());
        iEvTest->SetUNIXTime(evTest->GetUNIXTime());
        iEvTest->SetTBTime(evTest->GetTBTime());

        iEvTest->SetEventType(0);
        iEvTest->SetROIPixelID(evTest->GetROIPixelID());
        iEvTest->SetROIMusicID(evTest->GetROIMusicID());
        iEvTest->SetSignalValue(evTest->GetSignalValue());

        iEvTest->SetParametersFromTimestamp(timeAfterLaunch,data,data2);
        iEvTest->SetTelescopePointing(timeAfterLaunch,dataTilt);
        iEvTest->SetTelescopePointingRaw(timeAfterLaunch,dataTiltRaw);

        fOut->cd();
        tMTest->Fill();
    }

    fOut->Write();
    fOut->Close();
    fIn->Close();


    return 0;
}
