R__LOAD_LIBRARY(libExACT.so)
#include<TH1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <typeinfo>



int iLastPixHLED = -1;
int EventCounter = 0;
TLatex *text = 0;
TTree *tree = 0;
Event *ev;
long triggerTime;

TCanvas *c_disp = 0;
TH1I *hPixelTraceLED = 0;
BiFocal *bf = new BiFocal();

int MaxNofChannels = 256;
int LengthofPedestal = 200;
int SignalStart = 230;
int SignalWidth = 20;

void LoadEvents(std::string filename, std::string treeString);
void SetBranches(Event *evD, long *timeTrig, BiFocal *bfF);
void PlotEvent();
void ShowInfoAtCursor(int x, int y);

void FindBin(int Pixel_ID,int *row, int *col)
{
    // Calculate the SIAB number (0 to 15)
    int SIAB_Number = Pixel_ID / 16;

    // Calculate the pixel number within the SIAB (0 to 15)
    int SIAB_Pixel_Number = Pixel_ID % 16;

    // Calculate the row and column within the SIAB (0 to 3 for both)
    int SIAB_Pixel_Row = SIAB_Pixel_Number / 4;
    int SIAB_Pixel_Col = SIAB_Pixel_Number % 4;

    // Calculate the overall row and column
    *row = SIAB_Number / 4 * 4 + SIAB_Pixel_Row;
    *col = SIAB_Number % 4 * 4 + SIAB_Pixel_Col;
}

void TCameraPlotEvents(std::string filename,std::string treeString)
{
	c_disp = new TCanvas("Display","CameraPlot",950,1000);
	c_disp->Divide(2,2);
	LoadEvents(filename, treeString);
	ev = new Event();
	SetBranches(ev,&triggerTime, bf);
	PlotEvent();
}

Bool_t HandleInput()
{
	TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
	while (1)
	{
        	timer.TurnOn();
        	TString input = Getline("Type 'q' to exit, <return> to go on:");
        	timer.TurnOff();
        	if (input=="q\n") return kFALSE;
        	if (input=="\n") return kTRUE;
	}
	return kFALSE;
}

int FindPixel(int row, int col)
{
    // Calculate the SIAB number (0 to 15)
    int SIAB_Row = row / 4;
    int SIAB_Col = col / 4;
    int SIAB_Number = SIAB_Row * 4 + SIAB_Col;

    // Calculate the pixel number within the SIAB (0 to 15)
    int SIAB_Pixel_Row = SIAB_Row % 4;
    int SIAB_Pixel_Col = SIAB_Col % 4;
    int SIAB_Pixel_Number = SIAB_Pixel_Row * 4 + SIAB_Pixel_Col;

    // Calculate the pixel number (0 to 511)
    int Pixel_Number = SIAB_Number * 16 + SIAB_Pixel_Number;
    return Pixel_Number;
}

void LoadEvents(string NameofFile, std::string treeString)
{
    NameofFile = NameofFile;
    cout << "Loading file: " << NameofFile << endl;
    TFile *fO = new TFile(NameofFile.c_str(), "READ");
    tree = (TTree*)fO->Get(treeString.c_str());
}

void PlotTrace(int iPix)
{
    if(hPixelTraceLED == 0)
    {
        hPixelTraceLED = new TH1I("hPixelTraceLED","Pixel Trace",500,-0.5,499.5);
        hPixelTraceLED->SetStats(0);
        hPixelTraceLED->GetXaxis()->SetTitle("ADC sample");
        hPixelTraceLED->GetYaxis()->SetTitle("ADC counts");
    }

    hPixelTraceLED->Reset();
    TString title;
    title.Form("Trace of Pixel %i", iPix);
    hPixelTraceLED->SetTitle(title);
    std::vector<Int_t> trace = std::vector<Int_t>(512);

    trace  = ev->GetSignalValue(iPix);
    for(int k = 0; k<512; k++){
        hPixelTraceLED->SetBinContent(k+1,trace[k]);
    }

    c_disp->cd(3);
    hPixelTraceLED->Draw();

    gPad->Modified();
    gPad->Update();
}

void PixelClicked()
{
    int event = gPad->GetEvent();
    TObject *o = gPad->GetSelected();
    if (!o) return;
    if (!(o->InheritsFrom("TH2")))
       return;
    TH2F *h = (TH2F*)o;
    int px = gPad->GetEventX();
    int py = gPad->GetEventY();
    Float_t xx = gPad->AbsPixeltoX(px);
    Float_t yy = gPad->AbsPixeltoY(py);
    Float_t x = 0.5+gPad->PadtoX(xx);
    Float_t y = 0.5+gPad->PadtoY(yy);
    int pix = FindPixel((int)x,(int)y);
    if(pix!=iLastPixHLED)
    {
        int iLastPix = -1;
        ShowInfoAtCursor((int) x, (int) y);
    }
   if (event == 11)
    {
        PlotTrace(pix);
    }
}

void ShowInfoAtCursor(int x, int y)
{
    int MUSIC_column = x/2;
    int MUSIC_row = y/4;
    int MUSIC_ID = MUSIC_column+MUSIC_row*16;
    int MUSIC_Channel = y%4+4*(x%2);
    int PixID = MUSIC_row*8*16+MUSIC_column*8+MUSIC_Channel;

    TString statusline;
    statusline.Form("Event #: %i  MUSIC_ID: %i    MUSIC_Channel: %i    Pixel ID: %i",EventCounter, MUSIC_ID, MUSIC_Channel, PixID);
    if(text!=0)
        text->Delete();
    TLatex T1;
    text = T1.DrawLatexNDC(0.25,0.95,statusline.Data());
    gPad->Modified();
    gPad->Update();
}

void DrawMUSICBoundaries()
{
    TBox *b = new TBox(-.5,-0.5,3.5,1.5);
    b->SetFillStyle(0);
    b->SetLineColor(kRed);
    b->Draw();

    for(int i=0; i<MaxNofChannels/8; i++)
    {
        TBox *bn = (TBox*)b->Clone();
        bn->SetY1((i%8)*2-0.5);
        bn->SetY2((i%8)*2+1.5);
        bn->SetX1((i/8)*4-0.5);
        bn->SetX2((i/8)*4+3.5);
        bn->Draw();
    }
}

void PlotEvent()
{
    c_disp->cd(1);
    gPad->AddExec("ev","PixelClicked()");
    TH2F *hcam = new TH2F("hcam","",16,-0.5,15.5,16,-0.5,15.5);
    hcam->SetStats(0);
    hcam->Draw("colz");
    DrawMUSICBoundaries();
    int nEntries = tree->GetEntries();
    std::cout << "Total Number of Events: " << nEntries << std::endl;

    TH1F *hChg = new TH1F("hChg","Charge Distribution",410, 0, 4096);
    hChg->SetStats(0);
    hChg->GetXaxis()->SetTitle("Charge [ADC]");
    hChg->GetYaxis()->SetTitle("No. of Pixels");

    while(1)
    {
        tree->GetEntry(EventCounter);
        cout << "Event# " << EventCounter <<" is displayed." << endl;
        hcam->Reset();
        hChg->Reset();
        Int_t nEventsHLED = 0;
        ExtractedData *extractedDataHLED[MaxNofChannels];
        for (int i = 0; i<MaxNofChannels; i++){
            extractedDataHLED[i] = 0;
            extractedDataHLED[i] = new ExtractedData();
         }

        Pulse *pulse;
        for(int j = 0; j<MaxNofChannels; j++){
                pulse = new Pulse(ev->GetSignalValue(j));

                extractedDataHLED[j]->SetAmplitude(pulse->GetAmplitude());
                int nx, ny;
                FindBin(j,&nx,&ny);
                hcam->SetBinContent(nx+1,ny+1,extractedDataHLED[j]->GetAmplitude());
                hChg->Fill(extractedDataHLED[j]->GetAmplitude());

             }

    uint64_t convtime = ((uint64_t)ev->GetTBTime() / 100000000ULL);
    //Assuming timeTrig is in the format "seconds since epoch with 10^-7 second resolution"
    time_t epochTime = static_cast<time_t>(convtime); // Convert to seconds
    int microsec = convtime % 1000000;

    // Convert to struct tm in UTC
    tm* utcTime = gmtime(&epochTime);

    // Format and display the UTC time
    std::cout << "Trigger Time (UTC): ";
    std::cout << std::put_time(utcTime, "%Y-%m-%d %H:%M:%S") << "." << std::setfill('0') << std::setw(6) << microsec << std::endl;
    // cout << "time " << ev->GetTBTime() <<" is displayed." << endl;


        c_disp->cd(1);
	c_disp->SetTheta(90);
	c_disp->SetPhi(90);
        hcam->SetMinimum();
        hcam->SetMaximum();
        c_disp->cd(1)->Modified();
        c_disp->cd(1)->Update();

        c_disp->cd(2);
        hChg->Draw();
        c_disp->cd(2)->Modified();
        c_disp->cd(2)->Update();


        EventCounter++;
        if(EventCounter>=nEntries){
        	EventCounter = 0;
        }
        if(!HandleInput())
            break;
    }
}

void SetBranches(Event *evD, long *timeTrig, BiFocal *bfF)
{
    tree->SetBranchAddress("Events", &ev);
}

