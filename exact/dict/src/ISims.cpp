#include "ISims.h"


ISims::ISims(){

}
ISims::~ISims(){

}
double ISims::GeneralizedPoisson(double *x, double *par){
	float mu = par[0];
	float xTalkProb = par[1];
	float lambda = -1.0*TMath::Log(1.0-xTalkProb);

	float prob = mu*TMath::Power(mu+lambda*x[0],(x[0]-1.0))*TMath::Exp(-1.0*(mu+x[0]*lambda))/(TMath::Factorial((int)x[0]));

	return prob;
}

std::vector<std::vector<double>> ISims::OpticalXTalkExtraction(std::string filename){
	
	TGraph *gr = new TGraphErrors(filename.c_str(),"%lg %lg %lg %lg"," ");

	double *x,*y;

	x = gr->GetX();
	y = gr->GetEX();

	int nPoints = gr->GetN();
	std::vector<std::vector<double>> optXTalk = std::vector<std::vector<double>> (2,std::vector<double>(nPoints));

	for(int i = 0; i < nPoints; i++){
		optXTalk[0][i] = x[i];
		optXTalk[1][i] = y[i];
	}

	delete gr;

	return optXTalk;

}

std::vector<std::vector<int>> ISims::InvertCARETrace(std::vector<std::vector<int> *> traces){
	std::vector<std::vector<int>> invertedTrace = std::vector<std::vector<int>>(512, std::vector<int>(traces[0]->size()));
	for(int i = 0; i<traces.size();i++){
		for(int j = 0; j<traces[0]->size();j++){
			invertedTrace[i][j] = 4000 + (int)traces[i]->at(j)*-1.0;
		}
	}

	return invertedTrace;
}

std::vector<std::vector<int>> ISims::InvertCARETrace(std::vector<std::vector<int>> traces){
	std::vector<std::vector<int>> invertedTrace = std::vector<std::vector<int>>(traces.size(), std::vector<int>(traces[0].size()));
	for(int i = 0; i<traces.size();i++){
		for(int j = 0; j<traces[0].size();j++){
			invertedTrace[i][j] = 4000 + (int)traces[i][j]*-1.0;
		}
	}

	return invertedTrace;
}