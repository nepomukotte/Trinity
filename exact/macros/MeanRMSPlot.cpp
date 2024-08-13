R__LOAD_LIBRARY(libExACT.so)
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>

int iLastPixHLED = -1;
int EventCounter = 0;
TLatex *text = 0;
TTree *t = 0;
IEvent *ev;

BiFocal *bf = new BiFocal();
long triggerTime;

TCanvas *c_disp = 0;
int MaxNofChannels = 512;
int LengthofPedestal = 200;
Float_t SingleEventPedestal[512];
int SignalStart = 230;
int SignalWidth = 20;

void LoadEvents(std::string filename, std::string treeString);
void ShowInfoAtCursor(int x, int y);

// For plotting camera focal plane
void FindBin(int pixelID,int *nx, int *ny)
{
    int MUSIC_ID = pixelID/8;
    int MUSIC_Channel = pixelID%8; 
    int MUSIC_row = MUSIC_ID/16;
    int MUSIC_column = MUSIC_ID%16;

    *ny = 4*MUSIC_row+MUSIC_Channel%4; 
    *nx = 2*MUSIC_column+MUSIC_Channel/4; 
}

int FindPixel(int nx, int ny)
{
    int pixelID;
    pixelID = (nx*4+(ny%4)+(ny/4)*128);
    return pixelID;
}
void ShowInfoAtCursor(int x, int y)
{
    int MUSIC_column = x/2;
    int MUSIC_row = y/4;
    int MUSIC_ID = MUSIC_column+MUSIC_row*16;
    int MUSIC_Channel = y%4+4*(x%2);
    int PixID = MUSIC_row*8*16+MUSIC_column*8+MUSIC_Channel;

    TString statusline;
    statusline.Form("MUSIC_ID: %i      MUSIC_Channel: %i      Pixel ID: %i", MUSIC_ID, MUSIC_Channel, PixID);
    if(text!=0)
        text->Delete();
    TLatex T1;
    text = T1.DrawLatexNDC(0.25,0.95,statusline.Data());
    gPad->Modified();
    gPad->Update();
}
void DrawMUSICBoundaries()
{
    TBox *b = new TBox(-.5,-0.5,1.5,3.5);
    b->SetFillStyle(0);
    b->SetLineColor(kRed);
    b->Draw();
    //TPoint p;
    for(int i=0; i<MaxNofChannels/8; i++)
    {
        TBox *bn = (TBox*)b->Clone();
        bn->SetX1((i%16)*2-0.5);
        bn->SetX2((i%16)*2+1.5);
        bn->SetY1((i/16)*4-0.5);
        bn->SetY2((i/16)*4+3.5);
        bn->Draw();
    }
}
// grab event file and tree from terminal input
void LoadEvents(string NameofFile, std::string treeString)
{
    NameofFile = NameofFile;
    cout << "Loading file: " << NameofFile << endl;
    TFile *fO = new TFile(NameofFile.c_str(), "READ");
    t = (TTree*)fO->Get(treeString.c_str());
}

// main function: plots mean pedestal RMS
int MeanRMSPlot(std::string filename,std::string treeString){


	// where data for all pixels is stored
	ExtractedData *exD[512];
	
	// set up display and histograms
	c_disp = new TCanvas("Display","CameraPlot",950,1000);
	c_disp->Divide(1,2);
	LoadEvents(filename, "BiFocal");
	ev = new IEvent();
	t->SetBranchAddress("Events", &ev);


	
   c_disp->cd(1);
   TH2F *hcam = new TH2F("hcam","Average RMS over all events, May14th",32,-0.5,31.5,16,-0.5,15.5);
   hcam->SetStats(0);
   hcam->Draw("colz");
   DrawMUSICBoundaries();
   int nEntries = t->GetEntries();
   std::cout << "Total Number of Events: " << nEntries << std::endl;

	TH1 *h = new TH1F("h","RMS over all pixels and events",160,0,160);
	h->SetStats(0);
	h->GetXaxis()->SetTitle("RMS");
    h->GetYaxis()->SetTitle("No. of pixels*events");
	

std::vector<Double_t> accumulatedMean(MaxNofChannels, 0.0);

	for(int EventCounter = 0; EventCounter < nEntries; EventCounter++){

        t->GetEntry(EventCounter);
        cout << "Event# " << EventCounter <<" is displayed." << endl;
	//hcam->Reset();
   // h->Reset();
		Int_t nEventsHLED = 0;
        ExtractedData *extractedDataHLED[MaxNofChannels];
        for (int i = 0; i<MaxNofChannels; i++){
            extractedDataHLED[i] = 0;
            extractedDataHLED[i] = new ExtractedData();
         }

        Pulse *pulse;
       
		// set up variables for calculation
		Double_t PedRMS = 0;
		Double_t mean = 0;
        for (int j = 0; j<MaxNofChannels; j++){
			pulse = new Pulse(ev->GetSignalValue(j));
			extractedDataHLED[j]->SetPedestalRMS(pulse->GetPedestalRMS());
      		 accumulatedMean[j] += extractedDataHLED[j]->GetPedestalRMS();

			h->Fill(extractedDataHLED[j]->GetPedestalRMS()); // fill historgram with pedestal numbers
		}
	}
// Calculate the average meanPedestalRMS for each j
for (int j = 0; j < MaxNofChannels; j++) {
    accumulatedMean[j] /= nEntries;
   	int nx, ny;
    FindBin(j,&nx,&ny);
    hcam->SetBinContent(nx+1,ny+1, accumulatedMean[j]);
}

	// display all histograms
	c_disp->cd(1);
   hcam->SetMinimum(0);
   c_disp->cd(1)->Modified();
   c_disp->cd(1)->Update();

   c_disp->cd(2);
   h->Draw();
   c_disp->cd(2)->Modified();
   c_disp->cd(2)->Update();
    int numEntriesh = h->GetEntries();
    std::cout << "Number of Entries: " << numEntriesh << std::endl;

	return 0;
}
