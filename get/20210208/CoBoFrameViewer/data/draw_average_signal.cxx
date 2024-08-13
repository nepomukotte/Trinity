#include <GDataFrame.h>
#include <GDataSample.h>

#include <RVersion.h> // ROOT_VERSION_CODE
#include <TBenchmark.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TEntryListArray.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>
#include <TTree.h>
#include <iostream>

/**
 * Usage:
 * gSystem->Load("libget-root.so");
 * .L draw_channel_average_signal.cxx+
 * draw_channel_samples("CoBo_AsAd1_2013-08-23T15:33:49.836_0000.root", 2, 67, 100.0)
 * fit_ideal_signal()
 */
void draw_channel_average_signal(const std::string & filename, Int_t asadIndex=0, Int_t chipIndex=0, Int_t channelIndex=0, Double_t samplingFreq_MHz=100.0, Double_t triggerDelay_us=8., const char* drawOption="")
{
	gROOT->SetStyle("Plain");
	gSystem->Load("libget-root.so");

	gBenchmark->Start("draw_channel_average_signal");

	TFile* file = new TFile(filename.c_str());
	TTree *tree = (TTree*) file->Get("tree");
	if (! tree) return;

	tree->Print();

	gROOT->cd(); //Make sure that the new objects are memory resident

	// Create hist
	TH1F* hist[512];
	for (Int_t i=0; i< 512; i++) {
		hist[i] = new TH1F(TString::Format("hBucket%d", i), "", 4096, 0, 4096);
	}

	// Select event frames with channel
	tree->SetEntryList(0);
	tree->Draw(">>frameList", TString::Format("fHeader.fAsadIdx == %d && fChannels.fChanIdx == %d && fChannels.fAgetIdx == %d", asadIndex, channelIndex, chipIndex), "entrylistarray");
	TEntryListArray *frameList = (TEntryListArray*) gDirectory->Get("frameList");
	tree->SetEntryList(frameList);

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
		std::cout << "Loading frame # " << frameList->GetEntry(i) << std::endl;
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
			hist[sample->fBuckIdx]->Fill(sample->fValue);
		}
	}

	Double_t usPerBucket = 1.;
	if (samplingFreq_MHz > 0)
		usPerBucket = 1.0/samplingFreq_MHz;

	TGraphErrors* graph = new TGraphErrors(512);
	graph->SetName("gAverageSignal");
	graph->SetTitle(TString::Format("%s, AGET %d, Channel %d", filename.c_str(), chipIndex, channelIndex));
	for (Int_t i=0; i < 512; ++i) {
		Double_t t = usPerBucket*i;
		if (samplingFreq_MHz > 0)
			t = triggerDelay_us - usPerBucket*(511-i);
		graph->SetPoint(i, t, hist[i]->GetMean());
		graph->SetPointError(i, 0.0, hist[i]->GetRMS());
		hist[i]->Delete();
	}

	// Close file
	tree->ResetBranchAddresses();
	file->Close();

	gBenchmark->Show("draw_channel_average_signal");
	gBenchmark->Reset();

	TCanvas* canvas = 0;
	if (TString(drawOption).Contains("SAME")) {
		canvas = (TCanvas*) gROOT->FindObject("cChannelAverageSignal");
	}
	if (!canvas) {
		canvas = new TCanvas("cChannelAverageSignal", "", 800, 600);
		canvas->SetTopMargin(1.25);
		canvas->SetLeftMargin(0.10);
	}
	canvas->cd();
	gStyle->SetOptStat(kTRUE);
	gStyle->SetOptFit(kTRUE);

	graph->Draw(); //"H");
	graph->GetXaxis()->SetTitle("Time bucket [-]");
	if (samplingFreq_MHz > 0)
		graph->GetXaxis()->SetTitle("t - t_{0} [#mus]");
	graph->GetXaxis()->CenterTitle();
	graph->GetYaxis()->SetTitle("ADC value [-]");
	graph->GetYaxis()->CenterTitle();
	graph->GetYaxis()->SetTitleOffset(1.25);
	canvas->Modified();
	canvas->Update();
}

void fit_ideal_signal(Double_t tau_ns = 1014.)
{
	TGraphErrors* graph = (TGraphErrors*) gROOT->FindObject("gAverageSignal");
	if (!graph) return;
	TCanvas* canvas = (TCanvas*) gROOT->FindObject("cChannelAverageSignal");
	if (!canvas) return;
	canvas->cd();

	Long64_t n = graph->GetN();

	// Range for fit
	Double_t timeLow = TMath::MinElement(n, graph->GetX()), timeUp = TMath::MaxElement(n, graph->GetX());
	//Double_t timeLow =

	TF1*  func = new TF1("fAGET", "[0] + [1]*exp(-3*(x-[2])/[3])*((x-[2])/[3])**3*sin((x-[2])/[3])", timeLow, timeUp);
	func->SetLineColor(kRed);
	func->SetParNames("offset", "K", "t0", "#tau");
	// Initial value for pedestal
	Double_t adcMin = TMath::MinElement(n, graph->GetY());
	// Initial value for amplitude
	Double_t adcMax = TMath::MaxElement(n, graph->GetY());
	Long64_t indexMax = TMath::LocMax(n, graph->GetY());
	Double_t timeMax = (graph->GetX())[indexMax];
	func->SetParameters(adcMin, (adcMax-adcMin)/0.04409, timeMax - 1.166*tau_ns*1e-3, tau_ns*1e-3);
	graph->Fit(func, "RS");
	func->Print();
	func->Draw("SAME");

	// Compute peaking time
	Double_t time100pc = func->GetMaximumX(timeLow, timeUp);
	Double_t adc5pc = 0.95*func->GetParameter("offset") + 0.05* func->Eval(time100pc);
	Double_t time5pc = func->GetX(adc5pc, timeLow, time100pc);
	Double_t peakingTime_ns = 1000.*(time100pc - time5pc);
	std::cout << "Peaking time: " << peakingTime_ns << " ns" << std::endl;
	TBox* box = new TBox(time5pc, adc5pc, time100pc, func->Eval(time100pc));
	box->SetLineColor(kBlue);
	box->SetFillColor(kBlue);
	box->SetFillStyle(3004);
	box->Draw();
	TText* text = new TText(time5pc, func->Eval(time100pc), TString::Format("%.2lf ns", peakingTime_ns));
	text->SetTextColor(kBlue);
	text->SetTextSize(0.03);
	text->Draw();

	// Draw pedestals
	box = new TBox(timeLow, 0., timeUp, func->GetParameter("offset"));
	box->SetLineColor(kGreen);
	box->SetFillColor(kGreen);
	box->SetFillStyle(3005);
	box->Draw();

	// Plot residuals
	TGraph* graphRes = new TGraph(n);
	graphRes->SetName("gFitResiduals");
	graphRes->SetTitle("Signal fit residuals");
	for (Int_t i=0; i<n; i++) {
		double time, value;
		graph->GetPoint(i, time, value);
		graphRes->SetPoint(i, time, value - func->Eval(time));
	}
	canvas = new TCanvas("cSignalFitResiduals", "", 800, 600);
	canvas->SetTopMargin(1.25);
	canvas->SetLeftMargin(0.10);
	canvas->cd();

	graphRes->Draw();
	graphRes->GetXaxis()->SetTitle(graph->GetXaxis()->GetTitle());
	graphRes->GetXaxis()->CenterTitle();
	graphRes->GetYaxis()->SetTitle("Residual ADC value [-]");
	graphRes->GetYaxis()->CenterTitle();
	graphRes->GetYaxis()->SetTitleOffset(1.25);
	canvas->Modified();
	canvas->Update();
}
