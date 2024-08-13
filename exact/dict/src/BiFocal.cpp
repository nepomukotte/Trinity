#include "BiFocal.h"

BiFocal::BiFocal(){
	pxIDs = vector<unsigned short>(36);
	musicIDs = vector<unsigned short>(2);
}

BiFocal::~BiFocal(){
	
}

void BiFocal::SetROIPixelIDs(vector<int> rcvPxID){
	pxIDs.resize(rcvPxID.size());
	pxIDs = vector<unsigned short>(begin(rcvPxID), end(rcvPxID));
}

void BiFocal::SetTrigMUSICIDs(vector<int> rcvMusicID){
	musicIDs = vector<unsigned short>(begin(rcvMusicID),end(rcvMusicID));
}

void BiFocal::Copy(BiFocal *bf){
	pxIDs = bf->GetROIPixelIDs();
	musicIDs = bf->GetTrigMUSICIDs();
}

vector<unsigned short> BiFocal::GetROIPixelIDs(){
	return pxIDs;
}

vector<unsigned short> BiFocal::GetTrigMUSICIDs(){
	return musicIDs;
}
