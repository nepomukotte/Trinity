#ifndef ISIPM_H
#define ISIPM_H

#include <iostream>
#include <vector>

class ISiPM{
	public:
		ISiPM();
		ISiPM(int nPixels);
		~ISiPM();

		void SetGain(std::vector<double> values);
		void SetOptXTalk(std::vector<double> values);
		void SetTCorrection(std::vector<double> values);
		void SetAmplToPE(std::vector<double> values);
		void SetChgToPE(std::vector<double> values);

		std::vector<double> GetGain();
		std::vector<double> GetOptXTalk();
		std::vector<double> GetTCorrection();
		std::vector<double> GetAmplToPE();
		std::vector<double> GetChgToPE();

	private:
		std::vector <double> gain;
		std::vector <double> optXTalk;
		std::vector <double> gainTCorrection;
		std::vector <double> amplToPE;
		std::vector <double> chgToPE;
};

#endif // ISIPM_H