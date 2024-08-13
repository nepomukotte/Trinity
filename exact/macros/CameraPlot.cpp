R__LOAD_LIBRARY(libExACT.so)
#include<TH1.h>
#include <TTree.h>
#include <TCanvas.h>


int iLastPixHLED = -1;
TLatex *text = 0;
TTree *tree = 0;
TTree *tree_c = 0;
ExtractedData *ex[512];

vector<float> *correctionFactors = 0;//vector<float>(512,1.0);


BiFocal *bf = new BiFocal();
long triggerTime;

TCanvas *c_disp = 0;

int MaxNofChannels = 512;

void LoadEvents(std::string filename, std::string treeString);
void SetBranches(ExtractedData **exD, long *timeTrig, BiFocal *bfF);
void PlotEvent();
void ShowInfoAtCursor(int x, int y);

void FindBin(int pixelID,int *nx, int *ny)
{
    int MUSIC_ID = pixelID/8;
    int MUSIC_Channel = pixelID%8; 
    int MUSIC_row = MUSIC_ID/16;
    int MUSIC_column = MUSIC_ID%16;

    *ny = 4*MUSIC_row+MUSIC_Channel%4; 
    *nx = 2*MUSIC_column+MUSIC_Channel/4; 
}

void CameraPlot(std::string filename,std::string treeString){
	c_disp = new TCanvas("Display","CameraPlot",950,1000);
	c_disp->Divide(1,2);

	for(int i = 0; i<512; i++){
		ex[i] = new ExtractedData();
	}

	LoadEvents(filename, treeString);
	SetBranches(ex,&triggerTime, bf);
	PlotEvent();

}

Bool_t HandleInput(){
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
int FindPixel(int nx, int ny)
{
    int pixelID;
    pixelID = (nx*4+(ny%4)+(ny/4)*128);
    return pixelID;
}
void LoadEvents(string NameofFile, std::string treeString)
{
    NameofFile = NameofFile;
    cout << "Loading file: " << NameofFile << endl;
    TFile *fO = new TFile(NameofFile.c_str(), "READ");
    tree = (TTree*)fO->Get(treeString.c_str());
    tree_c = (TTree*)fO->Get("CorrectionFactor");
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
void PlotEvent()
{
    int EventCounter = 0;

    c_disp->cd(1);
    gPad->AddExec("ex","PixelClicked()");
    TH2F *hcam = new TH2F("hcam","",32,-0.5,31.5,16,-0.5,15.5);
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
    	//cout<<"Here"<<endl;
        tree->GetEntry(EventCounter);
        //tree_c->GetEntry(EventCounter);
        cout << "Event# " << EventCounter <<" is displayed." << endl;
        hcam->Reset();
        hChg->Reset();


        for(int i=0; i<512; i++)
        {
            int nx, ny;
            //Float_t tmp = 0;
            FindBin(i,&nx,&ny);
            /*for(int j=SignalStart; j<SignalStart+SignalWidth; j++)
            {
                tmp += Hled_struct.SignalValue[i][j] - SingleEventPedestal[i];
            }*/
            	//cout<<"Correction Factor "<<correctionFactors->at(i)<<endl;
              	hcam->SetBinContent(nx+1,ny+1,ex[i]->GetAmplitude());
				if(i==195){
					cout<<ex[195]->GetAmplitude()<<endl; 
				}
				     
           	hChg->Fill(ex[i]->GetAmplitude());
        }

        c_disp->cd(1);
        hcam->SetMinimum(0);
        //hcam->SetMaximum(2500);
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
void SetBranches(ExtractedData **exD, long *timeTrig, BiFocal *bfF)
{
    //tree->SetBranchAddress("TriggerTime", &timeTrig);
    for(int i = 0; i < 512; i++){
    	tree->SetBranchAddress(TString::Format("Pixel[%d]",i), &exD[i]);
    }
    //tree_c->SetBranchAddress("CorectionFactors", &correctionFactors);
    //tree->SetBranchAddress("TriggerRegion", bfF);
    //tHledEvents->SetBranchAddress("SignalValue", Hled_struct->SignalValue);
}

