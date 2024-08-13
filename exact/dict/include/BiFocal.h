#ifndef BIFOCAL_H
#define BIFOCAL_H

#include <TROOT.h>
#include "ReadConfiguration.h"
#include <vector>
#include "ExtractedData.h"

using namespace std;

class BiFocal{
public:
	BiFocal();
	~BiFocal();
	vector<unsigned short> GetROIPixelIDs();
	vector<unsigned short> GetTrigMUSICIDs();
	void SetROIPixelIDs(vector<int> rcvPxID);
	void SetTrigMUSICIDs(vector<int> rcvMusicID);
	void Copy(BiFocal *bf);

private:
	vector<unsigned short> pxIDs;
	vector<unsigned short> musicIDs;
};
#endif