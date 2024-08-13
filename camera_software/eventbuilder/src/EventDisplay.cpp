#include "EventDisplay.h"

using namespace std;

TTree *tHledEvents;
STRUCT_Hled Hled_struct;

Float_t Pedestal[MaxNofChannels];
Float_t SingleEventPedestal[MaxNofChannels];

int iLastPixHLED = -1;
TLatex *text = 0;
TCanvas *display1  = 0;
TH1I *hPixelTraceLED = 0;

void LoadSPB2Events(string NameofFile);
void SetHledBranch(STRUCT_Hled *Hled_struct);
void CalcHLEDPedestal();
void PlotSingleHLED();

void EventDisplay(std::string FileName)
{    
    display1 = new TCanvas("Display 3","Event Display",2000,2000);
    display1->Divide(2,2);

    LoadSPB2Events(FileName);
    SetHledBranch(&Hled_struct);    
    CalcHLEDPedestal();
    PlotSingleHLED();
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

int FindPixel(int nx, int ny)
{
    int pixelID;
    pixelID = (nx*4+(ny%4)+(ny/4)*128);
    return pixelID;
}

void FindBin(int pixelID,int *nx, int *ny)
{
    int MUSIC_ID = pixelID/8;
    int MUSIC_Channel = pixelID%8; 
    int MUSIC_row = MUSIC_ID/16;
    int MUSIC_column = MUSIC_ID%16;

    *ny = 4*MUSIC_row+MUSIC_Channel%4; 
    *nx = 2*MUSIC_column+MUSIC_Channel/4; 
}

void SetHledBranch(STRUCT_Hled *Hled_struct)
{
    tHledEvents->SetBranchAddress("TBTime", &Hled_struct->tb_time);
    tHledEvents->SetBranchAddress("CoBoTime", &Hled_struct->cobo_time);
    tHledEvents->SetBranchAddress("UNIXTime", &Hled_struct->unix_time);
    tHledEvents->SetBranchAddress("SignalValue", Hled_struct->SignalValue);
}

void LoadSPB2Events(string NameofFile)
{
    NameofFile = NameofFile + ".root";
    cout << "Loading file: " << NameofFile << endl;
    TFile *fO = new TFile(NameofFile.c_str(), "READ");
    tHledEvents = (TTree*)fO->Get("Hled_tree");
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
    //int iPixIndex = -1;
    for(int j=0; j<MaxTimeBucket; j++)
    {
        hPixelTraceLED->SetBinContent(j+1,Hled_struct.SignalValue[iPix][j]);
    }

    display1->cd(3);
    hPixelTraceLED->Draw();

    UInt_t min_val = Hled_struct.SignalValue[iPix][0];
    UInt_t max_val = Hled_struct.SignalValue[iPix][0];
    for(int i=0; i<MaxNofChannels; i++)
    {
    	if(min_val > Hled_struct.SignalValue[iPix][i]) {min_val = Hled_struct.SignalValue[iPix][i];}
    	if(max_val < Hled_struct.SignalValue[iPix][i]) {max_val = Hled_struct.SignalValue[iPix][i];}
    }

    TLine *lb = new TLine(SignalStart, min_val, SignalStart, max_val);
    lb->SetLineColor(kRed);
    lb->Draw();

    TLine *ub = new TLine(SignalStart+SignalWidth, min_val, SignalStart+SignalWidth, max_val);
    ub->SetLineColor(kRed);
    ub->Draw();

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
        PlotTrace(pix);
}

void CalcHLEDPedestal ()
{
    display1->cd(4);
    TH2F *hledDisplay4 = new TH2F("hledDisplay4","Pedestal Distribution",32,-0.5,31.5,16,-0.5,15.5);
    hledDisplay4->SetStats(0);
    hledDisplay4->Draw("colz");
    DrawMUSICBoundaries();

    int NofEvents = tHledEvents->GetEntries();
    for (int i=0; i<NofEvents; i++)
    {
        tHledEvents->GetEntry(i);
        for (int j=0; j<MaxNofChannels; j++)
        {
            Float_t tmp = 0;
            for(int k=0; k<LengthofPedestal; k++)
            {
                tmp += (1.0*Hled_struct.SignalValue[j][k]);
            }
            Pedestal[j] += (tmp/LengthofPedestal);
        }
    }
    for (int m=0; m<MaxNofChannels; m++)
    {
        int nx, ny;
        FindBin(m,&nx,&ny);
        Pedestal[m] = (Pedestal[m]/NofEvents);
        hledDisplay4->SetBinContent(nx+1, ny+1, Pedestal[m]);
        //hledDisplay4->SetMinimum(0);
        display1->cd(4)->Modified();
        display1->cd(4)->Update();
    }
};

void PlotSingleHLED()
{
    int EventCounter = 0;

    display1->cd(1);
    gPad->AddExec("ex","PixelClicked()");
    TH2F *hledDisplay1 = new TH2F("hledDisplay1","",32,-0.5,31.5,16,-0.5,15.5);
    hledDisplay1->SetStats(0);
    hledDisplay1->Draw("colz");
    DrawMUSICBoundaries();
    int totalnofentries = tHledEvents->GetEntries();
    std::cout << "Total Number of Events: " << totalnofentries << std::endl;

    TH1F *hledDisplay2 = new TH1F("hledDisplay2","Charge Distribution",200, 0, 20000);
    hledDisplay2->SetStats(0);
    hledDisplay2->GetXaxis()->SetTitle("Charge");
    hledDisplay2->GetYaxis()->SetTitle("Number of Pixels");

    while(1)
    {
        tHledEvents->GetEntry(EventCounter);
        cout << "Event# " << EventCounter <<" is displayed." << endl;
        hledDisplay1->Reset();
        hledDisplay2->Reset();

        for (int n=0; n<MaxNofChannels; n++)
        {
            Float_t tmp = 0;
            for(int m=0; m<LengthofPedestal; m++)
            {
                tmp += (1.0*Hled_struct.SignalValue[n][m]);
            }
            SingleEventPedestal[n] += (tmp/LengthofPedestal);
        }

        for(int i=0; i<MaxNofChannels; i++)
        {
            int nx, ny;
            Float_t tmp = 0;
            FindBin(i,&nx,&ny);
            for(int j=SignalStart; j<SignalStart+SignalWidth; j++)
            {
                tmp += Hled_struct.SignalValue[i][j] - SingleEventPedestal[i];
            }
            if(tmp > 0)
            {
            	hledDisplay1->SetBinContent(nx+1,ny+1,1);
            } else{
       	        hledDisplay1->SetBinContent(nx+1,ny+1,-1*tmp);
       	    }
            hledDisplay2->Fill(-1*tmp);
            SingleEventPedestal[i] = 0;
        }

        display1->cd(1);
        hledDisplay1->SetMinimum(0);
        display1->cd(1)->Modified();
        display1->cd(1)->Update();

        display1->cd(2);
        hledDisplay2->Draw();
        display1->cd(2)->Modified();
        display1->cd(2)->Update();

        EventCounter++;
        if(!HandleInput())
            break;
    }
}
