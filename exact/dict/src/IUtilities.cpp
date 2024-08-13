#include "IUtilities.h"

IUtilities::IUtilities(){

}
IUtilities::~IUtilities(){

}

uint64_t IUtilities::ConvertUTCtoUnix(const std::string& utcTimeStr) {
    // Convert UTC time string to tm structure
    std::tm utcTime = {};
    std::istringstream ss(utcTimeStr);
    ss >> std::get_time(&utcTime, "%Y-%m-%d %H:%M:%S");
    if (ss.fail()) {
        std::cerr << "Invalid UTC time format" << std::endl;
        return 0;
    }

    // Convert tm structure to Unix timestamp
    std::time_t unixTime = std::mktime(&utcTime);

    // Extract microseconds and convert to Unix timestamp format
    int microsec;
    ss >> microsec;

    // Calculate Unix timestamp with microseconds
    uint64_t unixTimestamp = static_cast<uint64_t>(unixTime) * 1000000ULL + microsec;

    return unixTimestamp;
}

int IUtilities::GetMaximumPixelID(std::vector<float>amplitude,int musicID){
	int iStartPixel, iEndPixel;

	if(musicID != -1){
		iStartPixel = musicID*8;
		iEndPixel = musicID*8+8;
	}else{
		iStartPixel = 0;
		iEndPixel = amplitude.size();
	}

	int max = 0;
	int maxID = 0;
	for(int i = iStartPixel; i<iEndPixel; i++){
		if(amplitude[i]>max){
			max = amplitude[i];
			maxID = i;
		}
	}

	return maxID;
}

int IUtilities::GetMaximumPixelID(std::vector<float>amplitude,std::vector<int> roiPixels){
	int iStartPixel, iEndPixel;


	int max = 0;
	int maxID = 0;

	if(roiPixels.size() == 0){
		iStartPixel = 0;
		iEndPixel = amplitude.size();

		for(int i = iStartPixel; i<iEndPixel; i++){
			if(amplitude[i]>max){
				max = amplitude[i];
				maxID = i;
			}
		}
	}else{
		for(int i = 0; i<roiPixels.size(); i++){
			if(amplitude[roiPixels[i]]>max){
				max = amplitude[roiPixels[i]];
				maxID = roiPixels[i];
			}
		}
	}
	return maxID;
}

float IUtilities::GetMaximum(std::vector<float>amplitude,int musicID){
	int iStartPixel, iEndPixel;

	if(musicID != -1){
		iStartPixel = musicID*8;
		iEndPixel = musicID*8+8;
	}else{
		iStartPixel = 0;
		iEndPixel = amplitude.size();
	}

	float max = 0;
	int maxID = 0;
	for(int i = iStartPixel; i<iEndPixel; i++){
		if(amplitude[i]>max){
			max = amplitude[i];
			maxID = i;
		}
	}

	return max;
}

float IUtilities::GetMaximum(std::vector<float>amplitude,std::vector<int> roiPixels){
	int iStartPixel, iEndPixel;


	float max = 0;
	int maxID = 0;

	if(roiPixels.size() == 0){
		iStartPixel = 0;
		iEndPixel = amplitude.size();

		for(int i = iStartPixel; i<iEndPixel; i++){
			if(amplitude[i]>max){
				max = amplitude[i];
				maxID = i;
			}
		}
	}else{
		for(int i = 0; i<roiPixels.size(); i++){
			if(amplitude[roiPixels[i]]>max){
				max = amplitude[roiPixels[i]];
				maxID = roiPixels[i];
			}
		}
	}
	return max;
}

int IUtilities::GetMinimumPixelID(std::vector<float>amplitude,int musicID){
	int iStartPixel, iEndPixel;

	if(musicID != -1){
		iStartPixel = musicID*8;
		iEndPixel = musicID*8+8;
	}else{
		iStartPixel = 0;
		iEndPixel = amplitude.size();
	}

	float min = 4096;
	int minID = 0;
	for(int i = iStartPixel; i<iEndPixel; i++){
		if(amplitude[i]<min){
			min = amplitude[i];
			minID = i;
		}
	}

	return minID;
}

int IUtilities::GetMinimumPixelID(std::vector<float>amplitude,std::vector<int> roiPixels){
	int iStartPixel, iEndPixel;


	float min = 4096;
	int minID = 0;

	if(roiPixels.size() == 0){
		iStartPixel = 0;
		iEndPixel = amplitude.size();

		for(int i = iStartPixel; i<iEndPixel; i++){
			if(amplitude[i]<min){
				min = amplitude[i];
				minID = i;
			}
		}
	}else{
		for(int i = 0; i<roiPixels.size(); i++){
			if(amplitude[roiPixels[i]]<min){
				min = amplitude[roiPixels[i]];
				minID = roiPixels[i];
			}
		}
	}
	return minID;
}

float IUtilities::GetMinimum(std::vector<float>amplitude,int musicID){
	int iStartPixel, iEndPixel;

	if(musicID != -1){
		iStartPixel = musicID*8;
		iEndPixel = musicID*8+8;
	}else{
		iStartPixel = 0;
		iEndPixel = amplitude.size();
	}

	float min = 4096;
	int minID = 0;
	for(int i = iStartPixel; i<iEndPixel; i++){
		if(amplitude[i]<min){
			min = amplitude[i];
			minID = i;
		}
	}

	return min;
}

float IUtilities::GetMinimum(std::vector<float>amplitude,std::vector<int> roiPixels){
	int iStartPixel, iEndPixel;


	float min = 4096;
	int minID = 0;

	if(roiPixels.size() == 0){
		iStartPixel = 0;
		iEndPixel = amplitude.size();

		for(int i = iStartPixel; i<iEndPixel; i++){
			if(amplitude[i]<min){
				min = amplitude[i];
				minID = i;
			}
		}
	}else{
		for(int i = 0; i<roiPixels.size(); i++){
			if(amplitude[roiPixels[i]]<min){
				min = amplitude[roiPixels[i]];
				minID = roiPixels[i];
			}
		}
	}
	return min;
}

float IUtilities::GetPedestalAverage(std::vector<int> trace, int tWindow,int pad, bool isUnCorrelated){
	std::vector<int> traceWindow = std::vector<int>(0);
	std::vector<int>::iterator itr;
	float average;
	float sum=0;
	int counter = 0;
	if(isUnCorrelated){
		for(int i = 0; i+tWindow<trace.size(); i+=(tWindow+pad)){
			traceWindow = std::vector<int>(trace.begin()+i,trace.begin()+i+tWindow);
			sum += std::accumulate(traceWindow.begin(),traceWindow.end(),0);
			counter++;
			traceWindow.clear();
		}
	}else{
		sum = std::accumulate(trace.begin(),trace.end(),0);
		counter = 1;
		tWindow = trace.size();
	}

	average = sum/(float)(counter*tWindow);
	return average;
}


std::vector<std::string> IUtilities::GetFilesInDirectory(std::string dirname, std::string extension){

	std::vector<std::string> files;
	std::string fileName;
	DIR *dir = opendir(dirname.c_str());

	if(dir){
		dirent *element;
		while ((element = readdir(dir))){
			fileName = element->d_name;
			if(fileName.length()>=extension.length()){
				if(fileName.compare(fileName.length() - extension.length(),extension.length(),extension) == 0){
					files.push_back(fileName);
				}
			}
		}

	}
	return files;
}

int IUtilities::GetHVChannel(int pixelID){
	int offset;
	offset = pixelID/256;
	while(pixelID>=128){
		pixelID = pixelID%128;

	}
	int hVChannel = (pixelID/32);

	return (hVChannel + offset*4);
}

double IUtilities::Interpolate(double x1, double y1, double x2, double y2, double xNew){

	double slope = (y2 - y1)/(x2 - x1);
    double offset = y2 - x2*slope;

   return slope*xNew + offset;
}
