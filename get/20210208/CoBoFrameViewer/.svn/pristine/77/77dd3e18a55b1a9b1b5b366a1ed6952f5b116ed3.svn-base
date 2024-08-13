#include <TBenchmark.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>

/**
 * Usage:
 * .L draw_multiplicity.cxx+
 * draw_multiplicity("CoBo_AsAd1_2013-08-23T15:33:49.836_0000.root")
 */
void draw_multiplicity(const std::string & filename)
{
	gROOT->SetStyle("Plain");


	gBenchmark->Start("draw_multiplicity");

	TFile* file = new TFile(filename.c_str());
	TTree *tree = (TTree*) file->Get("tree");
	if (! tree) return;

	//tree->Print();

	gROOT->cd(); //Make sure that the new objects are memory resident

	TCanvas* canvas = new TCanvas("cMultiplicity", "", 1);
	canvas->SetTopMargin(1.25);
	canvas->cd();
	gStyle->SetOptStat(kTRUE);
	gStyle->SetOptFit(kTRUE);

	// Fill hist
	Int_t frameCount = tree->Draw("Sum$(fHeader.fMult[]) >> hMultiplicity", "");
	std::cout << "Found " << frameCount << " frames(s)" << std::endl;

	// Close file
	file->Close();

	gBenchmark->Show("draw_multiplicity");
	gBenchmark->Reset();


	TH1D* hist = (TH1D*) gDirectory->Get("hMultiplicity");
	hist->GetXaxis()->SetTitle("AsAd board Multiplicity [-]");
	hist->GetXaxis()->CenterTitle();

	canvas->Modified();
	canvas->Update();
}

void draw_multiplicities(const std::string & filename)
{
	gROOT->SetStyle("Plain");

	gBenchmark->Start("draw_multiplicities");

	TFile* file = new TFile(filename.c_str());
	TTree *tree = (TTree*) file->Get("tree");
	if (! tree) return;

	//tree->Print();

	gROOT->cd(); //Make sure that the new objects are memory resident

	TCanvas* canvas = new TCanvas("cMultiplicities", "", 1024, 1024);
	canvas->cd();
	gStyle->SetOptStat(kTRUE);
	gStyle->SetOptFit(kTRUE);

	canvas->Divide(2, 2);

	// Loop over AGET chips
	for (Int_t agetIdx=0; agetIdx < 4; ++agetIdx) {
		TPad* pad = (TPad*) canvas->cd(1+agetIdx);
		pad->SetTopMargin(1.25);

		// Fill hist
		tree->Draw(TString::Format("fHeader.fMult[%d] >> hMultiplicity%d", agetIdx, agetIdx), "");
		TH1D* hist = (TH1D*) gDirectory->Get(TString::Format("hMultiplicity%d", agetIdx));
		hist->GetXaxis()->SetTitle(TString::Format("AGET %d Multiplicity [-]", agetIdx));
		hist->GetXaxis()->CenterTitle();

		pad->Modified();
	}

	// Close file
	file->Close();

	gBenchmark->Show("draw_multiplicities");
	gBenchmark->Reset();

	canvas->cd();
	canvas->Modified();
	canvas->Update();
}

void draw_windowOut(const std::string & filename)
{
	gROOT->SetStyle("Plain");


	gBenchmark->Start("draw_windowOut");

	TFile* file = new TFile(filename.c_str());
	TTree *tree = (TTree*) file->Get("tree");
	if (! tree) return;

	//tree->Print();

	gROOT->cd(); //Make sure that the new objects are memory resident

	TCanvas* canvas = new TCanvas("cWindowOut", "", 1);
	canvas->SetTopMargin(1.25);
	canvas->cd();
	gStyle->SetOptStat(kTRUE);
	gStyle->SetOptFit(kTRUE);

	// Fill hist
	Int_t frameCount = tree->Draw("fHeader.fWindowOut >> hWindowOut", "");
	std::cout << "Found " << frameCount << " frames(s)" << std::endl;

	// Close file
	file->Close();

	gBenchmark->Show("draw_windowOut");
	gBenchmark->Reset();


	TH1D* hist = (TH1D*) gDirectory->Get("hWindowOut");
	hist->GetXaxis()->SetTitle("Sliding Window Output [-]");
	hist->GetXaxis()->CenterTitle();

	canvas->Modified();
	canvas->Update();
}
