R__LOAD_LIBRARY(libExACT.so)
#include<TH1.h>
#include<TFitResult.h>


void FindBin(int pixelID,int *nx, int *ny)
{
    int MUSIC_ID = pixelID/8;
    int MUSIC_Channel = pixelID%8; 
    int MUSIC_row = MUSIC_ID/16;
    int MUSIC_column = MUSIC_ID%16;

    *ny = 4*MUSIC_row+MUSIC_Channel%4; 
    *nx = 2*MUSIC_column+MUSIC_Channel/4; 
}
            //Float_t tmp = 0;

int PlotHLEDDistribution(std::string fileName){
	ExtractedData *ex[512];
	vector<float> *correctionFactors = 0;
	unsigned long long timeTrig = 0;

	TFile *f = new TFile(fileName.c_str(),"READ");
	TTree *t = (TTree*)f->Get("HLED");
	//TTree *t_c = (TTree*)f->Get("CorrectionFactor");
	//t_c->SetBranchAddress("CorectionFactors",&correctionFactors);
	for(int i = 0; i < 512; i++){
		ex[i] = new ExtractedData();
		t->SetBranchAddress(TString::Format("Pixel[%d]",i),&ex[i]);
	}
	t->SetBranchAddress("TriggerTime",&timeTrig);

	int nx, ny;


	int nEntries = t->GetEntries();
	TH1 *h[512];

	/*for(int i = 0; i < 512; i++){
		h[i] = new TH1F(TString::Format("h[%d]",i),"",4096,0,4096);
	} */
	TH1 *h_col_L[5];
	TH1 *h_col_R[5];
	TH1 *h_all[5]; 
	TH1 *h_sum[5];
	TH1 *h_bP[5]; 

	TProfile2D *p_cam = new TProfile2D("p_cam","",32,-0.5,31.5,16,-0.5,15.5,0,4096);

	TProfile *h_time[5][4];

	for(int j = 0; j <5; j++){
		for(int i = 0; i<4; i++){
			h_time[j][i] = new TProfile(TString::Format("h_time[%d][%d]",j,i),"",120,0,120,0,4096);
		}
		
	}
	unsigned long long timeTrigStart;
	for(int i = 0; i<5;i++){
		h_col_L[i] = new TH1F(TString::Format("h_col_L[%d]",i),"",4096,0,4096);
		h_col_R[i] = new TH1F(TString::Format("h_col_R[%d]",i),"",4096,0,4096);
		h_all[i] = new TH1F(TString::Format("h_all[%d]",i),"",4096,0,4096);
		h_sum[i] = new TH1F(TString::Format("h_sum[%d]",i),"",4096,0,4096);
		h_bP[i] = new TH1F(TString::Format("h_bP[%d]",i),"",4096,0,4096);
	}
	for(int i = 0; i<nEntries; i++){
		t->GetEntry(i);
		//t_c->GetEntry(i);
		//cout<<timeTrig<<endl;
		if(i==0){
			timeTrigStart = timeTrig;
		}
		//cout<<(timeTrig-timeTrigStart)/1e8<<endl;
		if(i%100 == 0){
			cout<<"Processed: "<<i<<endl;
		}
		for(int j = 0; j<512; j++){
			FindBin(j,&nx,&ny);
			//cout<<timeTrig-timeTrigStart<<endl;
			
			if(nx>=20 && nx<= 23){
				h_col_R[i%5]->Fill(ex[j]->GetAmplitude());
				h_time[i%5][0]->Fill((timeTrig-timeTrigStart)/1e8,ex[j]->GetAmplitude());
			}else if(nx>=8 && nx<= 11){
				h_col_L[i%5]->Fill(ex[j]->GetAmplitude());
				h_time[i%5][1]->Fill((timeTrig-timeTrigStart)/1e8,ex[j]->GetAmplitude());
			}else if((nx>=16 && nx<= 19)&&(ny>=8 && ny<=11)){
				h_bP[i%5]->Fill(ex[j]->GetAmplitude());
			}else{
				h_time[i%5][2]->Fill((timeTrig-timeTrigStart)/1e8,ex[j]->GetAmplitude());
				h_all[i%5]->Fill(ex[j]->GetAmplitude());
				/*if(ex[j]->GetAmplitude()<1200 && i%5 == 3){
					cout<<j<<endl;
				}*/
			}
			if(i%5 ==3){
				// /cout<<nx<<endl;
				p_cam->Fill(nx,ny,ex[j]->GetAmplitude());
			}
			h_sum[i%5]->Fill(ex[j]->GetAmplitude());
			h_time[i%5][3]->Fill((timeTrig-timeTrigStart)/1e8,ex[j]->GetAmplitude());
		}
		
	}
	/*for(int j = 0; j<512; j++){
		//FindBin(j,&nx,&ny);
		h_avg->Fill(h[j]->GetMean());

	}*/

	cout<<"Entry No:\tRegion\tMean\tStd.Dev"<<endl;

	TCanvas *c = new TCanvas("c","",800,500);

	double meanL,meanR,meanAll,meanSum, meanBp;
	double rmsL,rmsR,rmsAll,rmsSum,rmsBp;

	double meanRawL, meanRawR, meanRawAll, meanRawSum;
	double rmsRawL, rmsRawR, rmsRawAll, rmsRawSum;

	TFitResultPtr r;
	for(int i = 2; i <3; i++){
		if(i==0){
			h_sum[i]->Draw();
			h_sum[i]->GetYaxis()->SetRangeUser(0,1400);
			h_sum[i]->GetYaxis()->SetTitle("N");
			h_sum[i]->GetXaxis()->SetTitle("ADC");
		}else{
			h_sum[i]->Draw("SAME");
			h_sum[i]->GetYaxis()->SetTitle("N");
			h_sum[i]->GetXaxis()->SetTitle("ADC");
		}

		r = h_sum[i]->Fit("gaus","QNS");
		meanRawSum = h_sum[i]->GetMean();
		meanSum = r->Parameter(1);
		rmsSum	 = r->Parameter(2);
		h_sum[i]->SetLineColor(kBlack);

		h_col_L[i]->Draw("SAME");
		h_col_L[i]->SetLineColor(kRed);
		r = h_col_L[i]->Fit("gaus","QNS");
		meanL = r->Parameter(1);
		rmsL	 = r->Parameter(2);
		meanRawL = h_col_L[i]->GetMean();

		h_col_R[i]->Draw("SAME");
		h_col_R[i]->SetLineColor(kGreen);
		r = h_col_R[i]->Fit("gaus","QNS");
		meanR = r->Parameter(1);
		rmsR	 = r->Parameter(2);
		meanRawR = h_col_R[i]->GetMean();

		h_all[i]->Draw("Same");
		h_all[i]->SetLineColor(kBlue);
		r = h_all[i]->Fit("gaus","QNS");
		meanAll = r->Parameter(1);
		rmsAll	 = r->Parameter(2);
		meanRawAll = h_all[i]->GetMean();

		h_bP[i]->Draw("SAME");
		h_bP[i]->SetLineColor(kMagenta);
		r = h_bP[i]->Fit("gaus","QNS");
		meanBp = r->Parameter(1);
		rmsBp = r->Parameter(2);

		cout<<i<<"\t\t"<<"Left\t"<<meanL<<"\t"<<rmsL<<"\t"<<meanRawL<<endl;
		cout<<i<<"\t\t"<<"Right\t"<<meanR<<"\t"<<rmsR<<"\t"<<meanRawR<<endl;
		cout<<i<<"\t\t"<<"Camera\t"<<meanAll<<"\t"<<rmsAll<<"\t"<<meanRawAll<<endl;
		cout<<i<<"\t\t"<<"Full\t"<<meanSum<<"\t"<<rmsSum<<"\t"<<meanRawSum<<endl;
		cout<<i<<"\t\t"<<"Bright\t"<<meanBp<<"\t"<<rmsBp<<"\t"<<endl;
		cout<<"----------------------------------------------"<<endl;

	}

	/*TCanvas *c_ts = new TCanvas("c_ts","",800,500);
	h_time[3][0]->Draw();
	h_time[3][1]->Draw("SAME");
	h_time[3][2]->Draw("SAME");
	h_time[3][3]->Draw("SAME");


	TCanvas *c_cam = new TCanvas("c_cam","", 1000,500);
	p_cam->Draw("COLZ");*/

	
	return 0;
}