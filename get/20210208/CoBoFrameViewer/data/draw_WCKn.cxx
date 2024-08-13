
#include <TGraphErrors.h>

/*
gSystem->Load("libget-root.so");
.L draw_average_signal.cxx+
.L draw_WCKn.cxx
draw_WCKn("ckw100_wckn025_wckn2100_triggerDelay800.root", 25., 100., 8.0)
draw_WCKn("ckw25_wckn06.25_wckn225_triggerDelay2360.root", 6.25, 25., 23.6)
 */
void draw_WCKn(const char* filename="ckw100_wckn025_wckn2100_triggerDelay800.root", double WCKn0=25., double WCKn2=100., double triggerDelay_us=8.)
{
	gSystem->Load("libget-root.so");

// ckw100_wckn025_wckn2100_triggerDelay800.root
// ckw25_wckn06.25_wckn225_triggerDelay1000.root
// ckw25_wckn06.25_wckn225_triggerDelay2360.root
draw_channel_average_signal("ckw100_wckn025_wckn2100_triggerDelay800.root", 0, 0, 0, WCKn0, triggerDelay_us);
TGraphErrors* gAsAd0 = (TGraphErrors*) gAverageSignal->Clone("gAsAd0");
draw_channel_average_signal("ckw100_wckn025_wckn2100_triggerDelay800.root", 2, 0, 0, WCKn2, triggerDelay_us);
TGraphErrors* gAsAd2 = (TGraphErrors*) gAverageSignal->Clone("gAsAd2");

TString title(TString::Format("triggerDelay=%g #mus (AGET 0, Channel 0)", triggerDelay_us));
gStyle->SetTitleBorderSize(0);
TCanvas* canvas = new TCanvas("cWCKn", "", 1);
canvas->cd();

gAsAd0->Draw("ALP");
gAsAd0->SetTitle(title);
gAsAd2->Draw("LPSAME");

gAsAd2->SetLineColor(kRed);
gAsAd2->SetMarkerColor(kRed);
gAsAd0->SetLineColor(kBlue);
gAsAd0->SetMarkerColor(kBlue);

TLegend *legend = new TLegend(0.1,0.7,0.4,0.8);
legend->SetFillStyle(0);
legend->SetBorderSize(0);
legend->SetHeader("");
legend->AddEntry(gAsAd0, TString::Format("WCKn=  %g MHz (AsAd0)", WCKn0), "l");
legend->AddEntry(gAsAd2, TString::Format("WCKn=  %g MHz (AsAd2)", WCKn2), "l");
legend->Draw();

}
