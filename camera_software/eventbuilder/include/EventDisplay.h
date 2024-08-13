#ifndef Format_Header
#define Format_Header

#include <TMath.h>
#include <TROOT.h>
#include <TTree.h>
#include <TLatex.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TTimer.h>
#include <TBox.h>
#include <TFile.h>
#include <TLine.h>
#include <TH2.h>
#include <Getline.h>

#include <iostream>


#define MaxNofChannels          512
#define MaxTimeBucket           512
#define LengthofPedestal        200
#define SignalStart         	230
#define SignalWidth         	20

using namespace std;
struct STRUCT_Hled
{
    ULong64_t tb_time;
    ULong64_t cobo_time;
    ULong64_t unix_time;
    UInt_t SignalValue[MaxNofChannels][MaxTimeBucket];
};

#endif
