#ifndef IHEALTHTOOLS_H
#define IHEALTHTOOLS_H

#include "Pulse.h"
#include "Event.h"
#include "IEvent.h"
#include "IFile.h"
#include <TString.h>
#include <TROOT.h>

class IHealthTools{
	public:
		IHealthTools();
		~IHealthTools();

		/*
		Obtains the full camera values for the specified measurement
		Measurement options are Pedestal, Pedestal RMS, Amplitude, 
		Charge and Peaking Time.

		All of the methods below require the same arguments.

		IFile *file: Object defined in ExACT. Contains the merged file and the 4 data Trees (HLED, BiFocal, Test, Forced)

		int entry: Entry to be extracted

		unsigned long long lTimeStart: start of time range

		unsigned long long lTimeEnd: end of time range

		std::string tree: tree name to be accessed. It is not case sensitive
		Can be (hled,bifocal,test,forced)

		float voltage: Checks the bias voltage from Auxiliary data. Tolerance
		is 0.1V.

		unsigned long long eventTime: Trigger time for the entry accessed.

		
		returns

		std:vector<float> with the value from the traces. 

		If the event is outside of time range or the camera was biased at a different voltage all elements are set to -1 and eventTime is set to 0

		*/

		static std::vector<float> CameraPedestal(IFile *file, int entry,unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string tree, float voltage, unsigned long long *eventTime,int iTWStart=200, int iTWEnd=300,float fTiltLow=-100,float fTiltHigh=100);
		static std::vector<float> CameraPedestalRMS(IFile *file, int entry,unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string tree, float voltage, unsigned long long *eventTime,int iTWStart=200, int iTWEnd=300,float fTiltLow=-100,float fTiltHigh=100);

		static std::vector<float> CameraAmplitude(IFile *file,int entry,unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string tree, float voltage, unsigned long long *eventTime,int iTWStart=200, int iTWEnd=300,float fTiltLow=-100,float fTiltHigh=100);

		static std::vector<float> CameraAmplitude(IFile *file,int entry, std::string tree, int iTWStart=200, int iTWEnd=300);
		static std::vector<float> CameraAmplitude(std::vector<std::vector<int>>, int iTWStart=200, int iTWEnd=300);

		static std::vector<float> CameraCharge(IFile *file,int entry,unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string tree, float voltage, unsigned long long *eventTime,int iTWStart=200, int iTWEnd=300,float fTiltLow=-100,float fTiltHigh=100);
		static std::vector<float> CameraPeakingTime(IFile *file,int entry,unsigned long long lTimeStart, unsigned long long lTimeEnd, std::string tree, float voltage, unsigned long long *eventTime,int iTWStart=200, int iTWEnd=300,float fTiltLow=-100,float fTiltHigh=100);

		static bool IsBadEvent(std::vector<std::vector<int>> traces, int iTWStart, int iTWEnd, int nPixelsTolerance = 20, int nSamplesTolerance = 10);
		static bool IsBadEvent(IFile *file,int entry, std::string treeName,int iTWStart, int iTWEnd, int nPixelsTolerance = 20, int nSamplesTolerance = 10);



	private:
		static bool IsTimeInRange(unsigned long long lTime,unsigned long long lTimeStart, unsigned long long lTimeEnd);
		static bool IsTiltInRange(float ftilt,float fTiltLow,float fTiltHigh);
};

#endif