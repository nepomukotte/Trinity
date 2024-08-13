#include <GDataFrame.h>
#include <GDataSample.h>

#include <RVersion.h> // ROOT_VERSION_CODE
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TEntryListArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>

/**
 * Usage:
 * gSystem->Load("libget-root.so");
 * .L draw_channel_samples.cxx+
 * draw_channel_samples("CoBo_AsAd1_2013-08-23T15:33:49.836_0000.root", 2, 67)
 */
void draw_channel_samples(const std::string & filename, Int_t chipIndex=0, Int_t channelIndex=0)
{
	gROOT->SetStyle("Plain");


	gBenchmark->Start("draw_channel_samples");

	TFile* file = new TFile(filename.c_str());
	TTree *tree = (TTree*) file->Get("tree");
	if (! tree) return;

	//tree->Print();

	gROOT->cd(); //Make sure that the new objects are memory resident

	// Select event frames with channel
	tree->SetEntryList(0);
	tree->Draw(">>frameList", TString::Format("fChannels.fChanIdx == %d && fChannels.fAgetIdx == %d", channelIndex, chipIndex), "entrylistarray");
	TEntryListArray *frameList = (TEntryListArray*) gDirectory->Get("frameList");
	tree->SetEntryList(frameList);

	// Create hist
	TH1F* hist = new TH1F("hChannelSamples", TString::Format("%s, AGET %d, Channel %d", filename.c_str(), chipIndex, channelIndex), 4096, 0, 4096);

#if ROOT_VERSION_CODE > 393216
	using GET::GDataFrame;
	using GET::GDataChannel;
	using GET::GDataSample;
#else
	using get::GDataFrame;
	using get::GDataChannel;
	using get::GDataSample;
#endif

	GDataFrame *frame = 0;
	tree->SetBranchAddress("GDataFrame", &frame);

	// Loop over event frames
	Int_t nentries = frameList->GetN();
	std::cout << "Looping over " << nentries << " frame(s)" << std::endl;
	for (Int_t i=0; i < nentries; ++i)
	{
		// Load frame
		//std::cout << "Loading frame # " << frameList->GetEntry(i) << std::endl;
		tree->GetEntry(frameList->GetEntry(i));
		TEntryListArray* channelList = frameList->GetSubListForEntry(frameList->GetEntry(i));

		// Get channel
		if (channelList->GetN() <= 0) continue;
		GDataChannel* channel = (GDataChannel*) frame->GetChannels()->At(channelList->GetEntry(0));

		// Loop over samples
		Int_t nsamples = channel->fSamples.GetEntries();
		//std::cout << "Looping over " << nsamples << " sample(s)" << std::endl;
		for (Int_t s = 0; s < nsamples; ++s)
		{
			GDataSample* sample = (GDataSample*) channel->fSamples.At(s);
			hist->Fill(sample->fValue);
		}
	}

	// Close file
	file->Close();

	gBenchmark->Show("draw_channel_samples");
	gBenchmark->Reset();

	TCanvas* canvas = new TCanvas("cChannelSamples", TString::Format("AGET %d, Channel %d", chipIndex, channelIndex).Data(), 1);
	canvas->SetTopMargin(1.25);
	canvas->cd();
	gStyle->SetOptStat(kTRUE);
	gStyle->SetOptFit(kTRUE);

	hist->Draw(); //"H"); // "H");
	hist->GetXaxis()->SetTitle("ADC value [-]");
	hist->GetXaxis()->CenterTitle();
	gPad->SetLogy(kTRUE);
	canvas->Modified();
	canvas->Update();
	hist->Fit("gaus");

	
}
