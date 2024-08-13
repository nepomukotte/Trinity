#include <iostream>
#include "ICalibration.h"

void CalibExtractor(int pixelID)
{
    ICalibration calibration("calibration.csv");

    //int pixelID = 1; // Replace with the desired pixel ID

    double amplitude = calibration.GetAmplitude(pixelID);
    double errAmplitude = calibration.GetErrAmplitude(pixelID);
    double charge = calibration.GetCharge(pixelID);
    double errCharge = calibration.GetErrCharge(pixelID);
    double relativeGain = calibration.GetRelativeGain(pixelID, 0, 4.7);

    std::cout << "Pixel ID: " << pixelID << std::endl;
    std::cout << "Amplitude: " << amplitude << std::endl;
    std::cout << "Error in Amplitude: " << errAmplitude << std::endl;
    std::cout << "Charge: " << charge << std::endl;
    std::cout << "Error in Charge: " << errCharge << std::endl;
    std::cout << "relativeGain: " << relativeGain << std::endl;
}
