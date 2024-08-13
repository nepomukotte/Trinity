#include <TROOT.h>
#ifndef IFHBEVENT_H
#define IFHBEVENT_H

#include "Event.h"

using namespace std;


// DataType class for all CT Events //

/* 

Please follow the following naming convention rules

Variables:
	- All public variables should start with capital letter I
		e.g. 
		int InumberOfPhotons
	- Only camel casing should be used. No underscores for variable names
		e.g. thisIsCameCasing this_Is_not
	- Please do not make all variables public if there is no need for them
	- Access to variable values should exclusively through Get/Set Methods

Methods:
	- Constructor method
		- The simplest constructor will initialize all variables to 0;
		- It is possible to overload for passing parameters
	- Camel Casing should be observed in naming all methods
	- Only ONE return statement for non-void methods
	- Getter and Setter Methods should be named Get[variableName] Set[variableName]
		e.g. GetNumberOfPhotons()
		- All Gets should return the datatype of the variable. 
		- All Sets should be void
	- Memory Address methods will be named Copy[variableName]
		e.g. int* CopyNumberOfPhotons()
		- Should always return a pointer to the variable in question.
	- Should only be relevant to data format itself. For analysis procedures please 
	use the include and source folders NOT in dict/
*/

class IEvent : public Event{
	public:
		/*EEvent DataType includes the following information

		ULong64_t tb_time; Time as recorded by the TB for the event in 10s of ns
		ULong64_t cobo_time; Time as recorded by the CoBo for the event
		ULong64_t unix_time; Time as recorded by the computer for the event
		Int_t eventType; Deprecated 
	
		vector<vector<Int_t>> signalValue; All traces for the event. contains a vecotr of 512x512
		vector<Int_t> roiMusicID; List of triggered music chips. Exactly 2 IDs are given
		vector<Int_t> roiPixelID; Region of Interest for the event 36 pixelIDs stored */

		IEvent(); //Constructor
		~IEvent(); //Destructor

	protected:
		// Telescope parameters //

		//float Iazimuth;
		//float Ialtitude;

		//...//

		// Calibration parameters //

		float IcorrectionFactor;

		//...//

	private:
		/*std::vector<std::string> headers;
		std::vector<std::vector<std::string>> data;
		std::vector<std::string> headers2;
		std::vector<std::vector<std::string>> data2;*/


		/*int timestampIndex;
		int timestamp2Index;

		int closestTimestamp1;
		int closestTimestamp2;*/

		// Parameters from the first file
		float latitude;
		float longitude;
		int altitude;
		float sunAzimuth;
		float sunElevation;
		float moonAzimuth;
		float moonElevation;
		float horizon;
		float azimuth;

		// Parameters from the second file
		int TrigEvent;
		int TempFlag;
		int emon1;
		int emon2;
		std::vector<float> hv;
		std::vector<float> hvc;
		std::vector<float> ucTemp;
		int CpuTemp;
		int CoboTemp;
		int RadTemp;
		float lvpsVol;
		float pumpVol;
		int lvpsCur;
		int pumpCur;
		std::vector<int> siabMPWR;
		std::vector<int> hvSW;
		std::vector<float> sipmTemp;

		// Tilt Sensor information
		float tiltAngle;
		float tiltAngleRaw;


	public:
		// void LoadFromFile(const std::string& filename1, const std::string& filename2);
		void SetParametersFromTimestamp(int timestamp,std::vector<std::vector<std::string>> data,std::vector<std::vector<std::string>> data2);
		void SetTelescopePointing(int timestamp, std::vector<std::vector<std::string>> data);
		void SetTelescopePointingRaw(int timestamp, std::vector<std::vector<std::string>> data);
		int GetclosestTimestamp1() const;
		int GetclosestTimestamp2() const;
		float GetLatitude() const;
		float GetLongitude() const;
		int GetAltitude() const;
		float GetSunAzimuth() const;
		float GetSunElevation() const;
		float GetMoonAzimuth() const;
		float GetMoonElevation() const;
		float GetHorizon() const;
		float GetAzimuth() const;
		int GetTrigEvent() const;
		int GetTempFlag() const;
		int GetEmon1() const;
		int GetEmon2() const;
		const std::vector<float>& Gethv() const;
		const std::vector<float>& Gethvc() const;
		const std::vector<float>& GetUCTemp() const;
		int GetCpuTemp() const;
		int GetCoboTemp() const;
		int GetRadTemp() const;
		float GetlvpsVol() const;
		float GetpumpVol() const;
		int GetlvpsCur() const;
		int GetpumpCur() const;
		const std::vector<int>& GetsiabMPWR() const;
		const std::vector<int>& GethvSW() const;
		const std::vector<float>& GetSiPMTemp() const;

		float GetTiltAngle() const;
		float GetTiltAngleRaw() const;

		static int FindClosestTimestampIndex(std::vector<unsigned long long> timeArray, unsigned long long eventTime);

	private:
    	int findClosestTimestamp(const std::vector<std::vector<std::string>>& data, int timestamp) const;

};
#endif
