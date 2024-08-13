#ifndef ISIMULATE_H
#define ISIMULATE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <TGraphErrors.h>
#include <TMath.h>

#define BREAKDOWN_V 38.7 //At 0deg C
class ISims{
	public:
		ISims();
		~ISims();
		static double GeneralizedPoisson(double *x, double *par);
		static std::vector<std::vector<double>> OpticalXTalkExtraction(std::string filename);
		static std::vector<std::vector<int>> InvertCARETrace(std::vector<std::vector<int> *> traces);
		static std::vector<std::vector<int>> InvertCARETrace(std::vector<std::vector<int>> traces);

};

#endif