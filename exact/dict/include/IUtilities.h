#ifndef IUTILITIES_H
#define IUTILITIES_H

#include <iostream>
#include <ctime>
#include <iomanip>
#include <vector>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <dirent.h>

class IUtilities{
	public:
		IUtilities();
		~IUtilities();


		/*
		Functions for random utilities to use as necessary for analysis:

		ConvertUTCtoUnix() : converts time

		*/

		static uint64_t ConvertUTCtoUnix(const std::string& utcTimeStr);





		static std::vector<std::string> GetFilesInDirectory(std::string dirname, std::string extension=".root");


		static int GetMaximumPixelID(std::vector<float> amplitude, int musicID=-1);
		static int GetMaximumPixelID(std::vector<float> amplitude, std::vector<int> roiPixels=std::vector<int>(0));
		static float GetMaximum(std::vector<float> amplitude, int musicID=-1);
		static float GetMaximum(std::vector<float> amplitude, std::vector<int> roiPixels=std::vector<int>(0));

		static int GetMinimumPixelID(std::vector<float> amplitude, int musicID=-1);
		static int GetMinimumPixelID(std::vector<float> amplitude, std::vector<int> roiPixels=std::vector<int>(0));
		static float GetMinimum(std::vector<float> amplitude, int musicID=-1);
		static float GetMinimum(std::vector<float> amplitude, std::vector<int> roiPixels=std::vector<int>(0));

		static float GetPedestalAverage(std::vector<int> trace, int tWindow=20,int pad=10, bool isUnCorrelated=true);

		static int GetHVChannel(int pixelID); 

		static double Interpolate(double x1, double y1, double x2, double y2, double xNew);
	
};

#endif
