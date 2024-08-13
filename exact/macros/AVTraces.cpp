R__LOAD_LIBRARY(libExACT.so)
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLatex.h>

int EventCounter = 0;
TTree *tree = 0;
IEvent *ev;

IEvent event;


BiFocal *bf = new BiFocal();
long triggerTime;

TCanvas *c_disp = 0;
TH1I *hAvgTraceLED = nullptr;
int MaxNofChannels = 512;


void LoadEvents(std::string filename, std::string treeString);
void SetBranches(Event *evD, long *timeTrig, BiFocal *bfF);
void PlotEvent();

void AVTraces(std::string filename, std::string treeString)
{
    c_disp = new TCanvas("Display", "CameraPlot", 1900, 1000);
    c_disp->Divide(1, 1);

    LoadEvents(filename, treeString);
    ev = new IEvent();
    SetBranches(ev, &triggerTime, bf);
    PlotEvent();
}

void LoadEvents(std::string NameofFile, std::string treeString)
{
    NameofFile = NameofFile;
    cout << "Loading file: " << NameofFile << endl;
    TFile *fO = new TFile(NameofFile.c_str(), "READ");
    tree = (TTree *)fO->Get(treeString.c_str());
}

void PlotEvent()
{
    int nEntries = tree->GetEntries();
    std::cout << "Total Number of Events: " << nEntries << std::endl;

    // Create a vector to store the average traces
    std::vector<TH1I*> avgTraces;
        int numPloted = 0;


	// loop through all events
    for (int EventCounter = 0; EventCounter < nEntries; EventCounter++)
    {
        tree->GetEntry(EventCounter);
        cout << "Event# " << EventCounter << " is displayed." << endl;

        // Create a new histogram for the average trace
        TH1I* hAvgTraceLED = new TH1I(Form("hAvgTraceLED_%d", EventCounter), "Average Trace", 500, -0.5, 499.5);
        hAvgTraceLED->SetStats(0);
        hAvgTraceLED->GetXaxis()->SetTitle("ADC sample");
        hAvgTraceLED->GetYaxis()->SetTitle("ADC counts");

        // Set axis label sizes
        hAvgTraceLED->GetXaxis()->SetLabelSize(0.04);
        hAvgTraceLED->GetYaxis()->SetLabelSize(0.04);
        
        // Set axis title sizes
        hAvgTraceLED->GetXaxis()->SetTitleSize(0.05);
        hAvgTraceLED->GetYaxis()->SetTitleSize(0.05);
        
        // Set histogram title size
        hAvgTraceLED->SetTitleSize(0.06);

        // Accumulate the average trace for each event
        std::vector<double> avgTrace(512, 0.0);
        for (int i = 0; i < MaxNofChannels; i++)
        {
            std::vector<Int_t> trace = ev->GetSignalValue(i);
            for (int j = 0; j < 512; j++)
            {
                avgTrace[j] += trace[j];
            }
        }
vector<float> hv;
vector<float> tmp;
int tmpFl;
     for (int i = 0; i < 8; i++)
        {
           hv = ev->Gethv();
           cout << "HV: " << hv[i] << endl;

        }
     /* for(int i = 0; i < 32; i++)
      {
      	tmp = ev->GetCamTemp();
      	cout << "temp: " << tmp[i] << endl;

      }

      	tmpFl = ev->GetTempFlag();
      	cout << "tempFlag: " << tmpFl << endl;
*/
        cout << "HV: " << hv[1] << endl;
        if (hv[1] <= 43.5)
            {
                    cout << "HV, plotted : " << hv[1] << endl;}

        // Normalize the accumulated trace by dividing by the total number of channels
        for (int k = 0; k < 512; k++)
        {
            avgTrace[k] /= MaxNofChannels;


	//	if(tmpFl == 1)
		{
	//	if(tmp[5] >= 0)
          //  if (hv[1] <= 43.5)
    {
            hAvgTraceLED->SetBinContent(k + 1, avgTrace[k]);
                                    numPloted++;  // Increment the counter for numPloted

            }
        }}

       

        // Store the average trace histogram in the vector
        avgTraces.push_back(hAvgTraceLED);

	// display plots
    c_disp->cd(1);
    hAvgTraceLED->GetXaxis()->SetRangeUser(220, 320);
    hAvgTraceLED->GetYaxis()->SetRangeUser(2000, 4000);
    hAvgTraceLED->SetTitle("Average Trace");
    hAvgTraceLED->Draw("SAME");
    c_disp->Modified();
    c_disp->Update();

    }
        cout << "Number of times hAvgTraceLED->SetBinContent was called: " << numPloted/512 << endl;

}


void SetBranches(Event *evD, long *timeTrig, BiFocal *bfF)
{
    tree->SetBranchAddress("Events", &ev);
}
