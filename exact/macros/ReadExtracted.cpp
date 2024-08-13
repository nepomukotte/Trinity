R__LOAD_LIBRARY(libExACT.so)
#include<TH1.h>
#include <TTree.h>
#include <TCanvas.h>


int ReadExtracted (){
	std::string filename =  "Extraction.root";
	ExtractedData *ex = new ExtractedData();

	TFile *f = new TFile(filename.c_str(), "READ");
	TTree *t = (TTree*)f->Get("t1");

	TH1 *h = new TH1F("h","hist",4096,0,4096);

	t->SetBranchAddress("Pixel[100]",&ex);

	for(int i = 0; i<t->GetEntries();i++){
		t->GetEntry(i);
		h->Fill(ex->GetAmplitude());
	}

	TCanvas *c = new TCanvas("c","Canvas",800,500);
	h->Draw();

	return 0;

}
