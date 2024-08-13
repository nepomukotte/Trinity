#ifndef GOLDPLATED_H
#define GOLDPLATED_H

#include <TROOT.h>
#include "ReadConfiguration.h"
#include <vector>
#include "ExtractedData.h"

using namespace std;

class GoldPlated : public ExtractedData{
public:
	GoldPlated();
	GoldPlated(int traceSize);
	~GoldPlated();

	vector<unsigned short> GetTrace();
	//vector<vector<unsigned short>> GetAllTraces();
	void SetTrace (vector<int> trace, int tStart);
	void CopyFromExtractedData(ExtractedData *ex);

protected:
	vector<unsigned short> signalValue;
};
#endif