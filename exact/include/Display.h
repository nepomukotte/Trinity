#ifndef Display_Class
#define Display_Class

#include<TROOT.h>
#include<TCanvas.h>
#include<TH1.h>
#include<TProfile2D.h>
#include <ReadConfiguration.h>
#include <Getline.h>
#include <TTimer.h>
#include <TApplication.h>

class Display{
public:
	Display(int argc, char **argv,ReadConfiguration *readConfig);
	~Display();
	void ShowPlotHLEDCoefficients(TProfile* p_hled);
	void ShowPlotHLEDCoefficientsDistribution(TH1* h_led[512], int pixNo);
protected:
	TApplication *app;
	int color[5] = {1,2,3,4,6};
};
#endif
