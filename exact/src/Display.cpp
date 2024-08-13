#include "Display.h"


Display::Display(int argc, char **argv,ReadConfiguration *readConfig){

	TROOT root("DisplayEvts","Display Results");
	app = new TApplication("App",&argc,argv);
}

Display::~Display(){

}

void Display::ShowPlotHLEDCoefficients(TProfile *p_hled){
	TCanvas *c = new TCanvas("c_hled", "Canvas",800,500);
	//Only 5 amplitudes on HLED
	c->cd(1);
	/*for(int i = 0; i < 5; i++){
		if(i==0){
			p_hled[i]->Draw("HIST P");
		}
		else{
			p_hled[i]->Draw("HIST P Same");
		}
		p_hled[i]->SetMarkerColor(color[i]);
		p_hled[i]->SetMarkerStyle(3);
		p_hled[i]->SetMarkerSize(3);
		
	}*/
	p_hled->Draw("");
	p_hled->GetZaxis()->SetRangeUser(0,5);

	TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
	timer.TurnOn();
	TString input = Getline("Type <return> to go on: ");
	timer.TurnOff();
}

void Display::ShowPlotHLEDCoefficientsDistribution(TH1 *h_led[512], int pixNo){
	TCanvas *c = new TCanvas("c_hled", "Canvas",800,500);
	//Only 5 amplitudes on HLED
	c->cd(1);
	/*for(int i = 0; i < 5; i++){
		if(i==0){
			p_hled[i]->Draw("HIST P");
		}
		else{
			p_hled[i]->Draw("HIST P Same");
		}
		p_hled[i]->SetMarkerColor(color[i]);
		p_hled[i]->SetMarkerStyle(3);
		p_hled[i]->SetMarkerSize(3);
		
	}*/
	h_led[pixNo]->Draw("");
	//h_led->GetZaxis()->SetRangeUser(0,5);

	TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
	timer.TurnOn();
	TString input = Getline("Type <return> to go on: ");
	timer.TurnOff();
}
