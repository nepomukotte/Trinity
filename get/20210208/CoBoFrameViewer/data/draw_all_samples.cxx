#include <TBenchmark.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>

void draw_all_samples(const std::string & filename)
{
	gROOT->SetStyle("Plain");


	gBenchmark->Start("draw_all_samples");

	TFile* file = new TFile(filename.c_str());
	TTree *tree = (TTree*) file->Get("tree");
	if (! tree) return;

	//tree->Print();

	gROOT->cd(); //Make sure that the new objects are memory resident

	// Create hist
	TH1F* hist = new TH1F("hSamples", TString::Format("%s", filename.c_str()), 4096, 0, 4096);

	// Fill hist
	Int_t sampleCount = tree->Draw("fSamples.fValue >> hSamples", "", "goff");
	std::cout << "Found " << sampleCount << " sample(s)" << std::endl;

	// Close file
	file->Close();

	gBenchmark->Show("draw_all_samples");
	gBenchmark->Reset();

	TCanvas* canvas = new TCanvas("cSamples", "", 1);
	canvas->SetTopMargin(1.25);
	canvas->cd();
	gStyle->SetOptStat(kTRUE);
	gStyle->SetOptFit(kTRUE);

	hist->Draw(); //"H");
	hist->GetXaxis()->SetTitle("ADC value [-]");
	hist->GetXaxis()->CenterTitle();
	gPad->SetLogy(kTRUE);
	canvas->Modified();
	canvas->Update();
	hist->Fit("gaus");
}
