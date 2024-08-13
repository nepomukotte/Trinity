#include "IHealthTools.h"

IHealthTools::IHealthTools(){

}
IHealthTools::~IHealthTools(){

}

bool IHealthTools::IsTimeInRange(unsigned long long lTime,unsigned long long lTimeStart, unsigned long long lTimeEnd){

	return (lTime>=lTimeStart && lTime<=lTimeEnd);
}
bool IHealthTools::IsTiltInRange(float fTilt,float fTiltLow,float fTiltHigh){

	return (fTilt>=fTiltLow && fTilt<=fTiltHigh);
}

std::vector<float> IHealthTools::CameraPeakingTime(IFile *file, int entry,unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string treeName, float voltage, unsigned long long *eventTime,int iTWStart, int iTWEnd,float fTiltLow,float fTiltHigh){

	IEvent *event = 0;
	std::vector<float> hv ;
	float tiltAngle;
	unsigned long long timeStamp;
	int nPixelsCamera = 512;

	std::vector<float> pkTime = std::vector<float>(nPixelsCamera);
	
	//TH2 *h = new TH2D("h","",nBinsX,lTimeStart,lTimeEnd,nBinsY,0,nADC);

	Pulse *p;

	int selectedTree = IFile::GetTreeID(treeName);

	switch(selectedTree){
		case 0:
			file->treeHLED->SetBranchAddress("Events",&event);

			file->treeHLED->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			//cout<<"HV: "<<hv[0]<<endl;
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pkTime[j] = p->GetTimePeak();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pkTime[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 1:
			file->treeBiFocal->SetBranchAddress("Events",&event);

			file->treeBiFocal->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();

			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pkTime[j] = p->GetTimePeak();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pkTime[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 2:
			file->treeTest->SetBranchAddress("Events",&event);

			file->treeTest->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();

			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pkTime[j] = p->GetTimePeak();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pkTime[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 3:
			file->treeForced->SetBranchAddress("Events",&event);

			file->treeForced->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();

			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pkTime[j] = p->GetTimePeak();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pkTime[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		break;
	}
	delete event;
	//delete p;

	return pkTime;
}
std::vector<float> IHealthTools::CameraAmplitude(IFile *file, int entry,unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string treeName, float voltage, unsigned long long *eventTime,int iTWStart, int iTWEnd,float fTiltLow,float fTiltHigh){

	IEvent *event = 0;
	std::vector<float> hv ;
	unsigned long long timeStamp;
	float tiltAngle;
	int nPixelsCamera = 512;

	std::vector<float> amplitude = std::vector<float>(nPixelsCamera);
	
	//TH2 *h = new TH2D("h","",nBinsX,lTimeStart,lTimeEnd,nBinsY,0,nADC);

	Pulse *p;	

	int selectedTree = IFile::GetTreeID(treeName);

	switch(selectedTree){
		case 0:
			file->treeHLED->SetBranchAddress("Events",&event);

			file->treeHLED->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle() - event->GetHorizon();
			//cout<<"HV: "<<hv[0]<<endl;
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.2 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					amplitude[j] = p->GetAmplitude();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					amplitude[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 1:
			file->treeBiFocal->SetBranchAddress("Events",&event);

			file->treeBiFocal->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle() - event->GetHorizon();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					amplitude[j] = p->GetAmplitude();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					amplitude[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 2:
			file->treeTest->SetBranchAddress("Events",&event);

			file->treeTest->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle() - event->GetHorizon();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					amplitude[j] = p->GetAmplitude();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					amplitude[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 3:
			file->treeForced->SetBranchAddress("Events",&event);

			file->treeForced->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle() - event->GetHorizon();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					amplitude[j] = p->GetAmplitude();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					amplitude[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		default:
			break;
	}
	
	delete event;
	return amplitude;
}

std::vector<float> IHealthTools::CameraAmplitude(IFile *file, int entry, std::string treeName, int iTWStart, int iTWEnd){

	Event *event = 0;
	std::vector<float> hv ;
	unsigned long long timeStamp;
	float tiltAngle;
	int nPixelsCamera = 512;

	std::vector<float> amplitude = std::vector<float>(nPixelsCamera);
	
	//TH2 *h = new TH2D("h","",nBinsX,lTimeStart,lTimeEnd,nBinsY,0,nADC);

	Pulse *p;	

	int selectedTree = IFile::GetTreeID(treeName);

	switch(selectedTree){
		case 0:
			file->treeHLED->SetBranchAddress("Events",&event);

			file->treeHLED->GetEntry(entry);
			
			//cout<<"HV: "<<hv[0]<<endl;
			
			for(int j= 0; j<nPixelsCamera; j++){
				p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
				amplitude[j] = p->GetAmplitude();
				delete p;
			}
			
			break;
		case 1:
			file->treeBiFocal->SetBranchAddress("Events",&event);

			file->treeBiFocal->GetEntry(entry);
			
			for(int j= 0; j<nPixelsCamera; j++){
				p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
				amplitude[j] = p->GetAmplitude();
				delete p;
			}
			break;
		case 2:
			file->treeTest->SetBranchAddress("Events",&event);

			file->treeTest->GetEntry(entry);
			
			for(int j= 0; j<nPixelsCamera; j++){
				p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
				amplitude[j] = p->GetAmplitude();
				delete p;
			}
			break;
		case 3:
			file->treeForced->SetBranchAddress("Events",&event);

			file->treeForced->GetEntry(entry);
			
			for(int j= 0; j<nPixelsCamera; j++){
				p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
				amplitude[j] = p->GetAmplitude();
				delete p;
			}
			break;
		default:
			break;
	}
	
	delete event;
	return amplitude;
}

std::vector<float> IHealthTools::CameraAmplitude(std::vector<std::vector<int>>traces, int iTWStart, int iTWEnd){

	
	int nPixelsCamera = 512;

	std::vector<float> amplitude = std::vector<float>(nPixelsCamera);
	
	//TH2 *h = new TH2D("h","",nBinsX,lTimeStart,lTimeEnd,nBinsY,0,nADC);

	Pulse *p;	

	for(int j= 0; j<nPixelsCamera; j++){
		p = new Pulse(traces[j],iTWStart,iTWEnd);
		amplitude[j] = p->GetAmplitude();
		delete p;
	}
			
	return amplitude;
}

std::vector<float> IHealthTools::CameraPedestal(IFile *file,int entry,unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string treeName, float voltage, unsigned long long *eventTime,int iTWStart, int iTWEnd,float fTiltLow,float fTiltHigh){

	IEvent *event = 0;
	std::vector<float> hv ;
	unsigned long long timeStamp;
	float tiltAngle;
	int nPixelsCamera = 512;

	std::vector<float> pedestal = std::vector<float>(nPixelsCamera);
	
	//TH2 *h = new TH2D("h","",nBinsX,lTimeStart,lTimeEnd,nBinsY,0,nADC);

	Pulse *p;

	int selectedTree = IFile::GetTreeID(treeName);

	switch(selectedTree){
		case 0:
			file->treeHLED->SetBranchAddress("Events",&event);

			file->treeHLED->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			//cout<<"Tilt Angle: "<<tiltAngle<<endl;
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pedestal[j] = (float)p->GetPedestal();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pedestal[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 1:
			file->treeBiFocal->SetBranchAddress("Events",&event);

			file->treeBiFocal->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pedestal[j] = p->GetPedestal();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pedestal[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 2:
			file->treeTest->SetBranchAddress("Events",&event);

			file->treeTest->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pedestal[j] = p->GetPedestal();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pedestal[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 3:
			file->treeForced->SetBranchAddress("Events",&event);

			file->treeForced->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pedestal[j] = p->GetPedestal();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pedestal[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		default:
			break;
	}
	
	delete event;
	//delete p;
	return pedestal;
}

std::vector<float> IHealthTools::CameraCharge(IFile *file,int entry,unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string treeName, float voltage, unsigned long long *eventTime,int iTWStart, int iTWEnd,float fTiltLow,float fTiltHigh){
	IEvent *event = 0;
	std::vector<float> hv ;
	unsigned long long timeStamp;
	int nPixelsCamera = 512;

	float tiltAngle;

	std::vector<float> charge = std::vector<float>(nPixelsCamera);
	
	//TH2 *h = new TH2D("h","",nBinsX,lTimeStart,lTimeEnd,nBinsY,0,nADC);

	Pulse *p;
	

	int selectedTree = IFile::GetTreeID(treeName);

	switch(selectedTree){
		case 0:
			file->treeHLED->SetBranchAddress("Events",&event);

			file->treeHLED->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			//cout<<"HV: "<<hv[0]<<endl;
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					charge[j] = p->GetCharge();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					charge[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 1:
			file->treeBiFocal->SetBranchAddress("Events",&event);

			file->treeBiFocal->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					charge[j] = p->GetCharge();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					charge[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 2:
			file->treeTest->SetBranchAddress("Events",&event);

			file->treeTest->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					charge[j] = p->GetCharge();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					charge[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 3:
			file->treeForced->SetBranchAddress("Events",&event);

			file->treeForced->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					charge[j] = p->GetCharge();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					charge[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		default:
			break;
	}
	
	delete event;
	return charge;
}


std::vector<float> IHealthTools::CameraPedestalRMS(IFile *file, int entry, unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string treeName, float voltage, unsigned long long *eventTime,int iTWStart, int iTWEnd,float fTiltLow,float fTiltHigh){
	
	IEvent *event = 0;
	std::vector<float> hv ;
	unsigned long long timeStamp;
	int nPixelsCamera = 512;
	float tiltAngle;

	std::vector<float> pedestalRMS = std::vector<float>(nPixelsCamera);
	
	//TH2 *h = new TH2D("h","",nBinsX,lTimeStart,lTimeEnd,nBinsY,0,nADC);
	Pulse *p;

	

	int selectedTree = IFile::GetTreeID(treeName);

	switch(selectedTree){
		case 0:
			file->treeHLED->SetBranchAddress("Events",&event);

			file->treeHLED->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			//cout<<"HV: "<<hv[0]<<endl;
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pedestalRMS[j] = p->GetPedestalRMS();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pedestalRMS[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 1:
			file->treeBiFocal->SetBranchAddress("Events",&event);

			file->treeBiFocal->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<=0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pedestalRMS[j] = p->GetPedestalRMS();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pedestalRMS[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 2:
			file->treeTest->SetBranchAddress("Events",&event);

			file->treeTest->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pedestalRMS[j] = p->GetPedestalRMS();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pedestalRMS[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		case 3:
			file->treeForced->SetBranchAddress("Events",&event);

			file->treeForced->GetEntry(entry);
			timeStamp = event->GetTBTime();
			hv = event->Gethv();
			tiltAngle = event->GetTiltAngle();
			
			if(timeStamp>=lTimeStart && timeStamp<=lTimeEnd && (TMath::Abs(hv[0]-voltage))<0.1 && IsTiltInRange(tiltAngle,fTiltLow,fTiltHigh)){
				*eventTime = timeStamp;
				for(int j= 0; j<nPixelsCamera; j++){
					p = new Pulse(event->GetSignalValue(j),iTWStart,iTWEnd);
					pedestalRMS[j] = p->GetPedestalRMS();
					delete p;
				}
			}else{
				for(int j= 0; j<nPixelsCamera; j++){
					pedestalRMS[j] = -1;
				}
				*eventTime = 0;
			}
			break;
		default:
			break;
	}
	delete event;
	//delete p;
	return pedestalRMS;
}

bool IHealthTools::IsBadEvent(std::vector<std::vector<int>> traces, int iTWStart, int iTWEnd, int nPixelsTolerance, int nSamplesTolerance){
	int nPixelsCamera = 512;

	int nTraceSamplesBad = 0;
	int nPixelsBad = 0;

	bool isBadEvent = false;
	int j = 0;
	for(int i = 0; i<nPixelsCamera; i++){
		j = iTWStart;
		while(j<iTWEnd && nTraceSamplesBad < nSamplesTolerance){
			if(traces[i][j] == traces[i][j+1]){
				nTraceSamplesBad++;
			}
			else{
				nTraceSamplesBad = 0;
			}
			j++;
		}
		if(nTraceSamplesBad == nSamplesTolerance){
			nPixelsBad++;
		}
		nTraceSamplesBad = 0;	
	}

	if(nPixelsBad > nPixelsTolerance){
		isBadEvent = true;
	}

	return isBadEvent;
}

bool IHealthTools::IsBadEvent(IFile *file,int entry, std::string treeName, int iTWStart, int iTWEnd, int nPixelsTolerance, int nSamplesTolerance){
	int nPixelsCamera = 512;

	int nTraceSamplesBad = 0;
	int nPixelsBad = 0;

	bool isBadEvent = false;
	int selectedTree = IFile::GetTreeID(treeName);
	TTree *treeIn;
	switch(selectedTree){
		case 0:
			treeIn = file->treeHLED;
			break;
		case 1:
			treeIn = file->treeBiFocal;
			break;
		case 2:
			treeIn = file->treeTest;
			break;
		case 3:
			treeIn = file->treeForced;
			break;
		default:
			cout<<"No such tree"<<endl;
			break;
	}
	int j = 0;
	IEvent *ev;
	std::vector<std::vector<int>> traces;
	ev = new IEvent();
	treeIn->SetBranchAddress("Events", &ev);
	treeIn->GetEntry(entry);

	traces = ev->GetSignalValue();

	for(int i = 0; i<nPixelsCamera; i++){
		j = iTWStart;
		while(j<iTWEnd && nTraceSamplesBad < nSamplesTolerance){
			if(traces[i][j] == traces[i][j+1]){
				nTraceSamplesBad++;
			}
			else{
				nTraceSamplesBad = 0;
			}
			j++;
		}
		if(nTraceSamplesBad == nSamplesTolerance){
			nPixelsBad++;
		}
		nTraceSamplesBad = 0;	
	}

	if(nPixelsBad > nPixelsTolerance){
		isBadEvent = true;
	}
	delete ev;
	//delete treeIn;
	return isBadEvent;
}

