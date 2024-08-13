#include <iostream>
#include "IEvent.h"

void AuxExtractor(std::string filename1, std::string filename2, int timestamp) {

    // Create an instance of the IEvent class
    IEvent event;

    // Load the files
    event.LoadFromFile(filename1, filename2);

    // Set the parameters based on the timestamp
    event.SetParametersFromTimestamp(timestamp);

    // Print the parameters
    std::cout << "\nClosest timestamp on 1st file: " << event.GetclosestTimestamp1() << std::endl;
    std::cout << "Latitude: " << event.GetLatitude() << std::endl;
    std::cout << "Longitude: " << event.GetLongitude() << std::endl;
    std::cout << "Altitude: " << event.GetAltitude() << std::endl;
    std::cout << "Sun Azimuth: " << event.GetSunAzimuth() << std::endl;
    std::cout << "Sun Elevation: " << event.GetSunElevation() << std::endl;
    std::cout << "Moon Azimuth: " << event.GetMoonAzimuth() << std::endl;
    std::cout << "Moon Elevation: " << event.GetMoonElevation() << std::endl;
    std::cout << "Horizon: " << event.GetHorizon() << std::endl;

    std::cout << "\nClosest timestamp on 2nd file: " << event.GetclosestTimestamp2() << std::endl;
    std::cout << "iTrigEvent: " << event.GetTrigEvent() << std::endl;
    std::cout << "iTempFlag: " << event.GetTempFlag() << std::endl;
    std::cout << "emon1: " << event.GetEmon1() << std::endl;
    std::cout << "emon2: " << event.GetEmon2() << std::endl;
    const std::vector<float>& hv = event.Gethv();
    std::cout << "Ihv: ";
    for (float value : hv) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    const std::vector<float>& hvc = event.Gethvc();
    std::cout << "hvc: ";
    for (float value : hvc) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    const std::vector<float>& camTemp = event.GetCamTemp();
    std::cout << "CamTemp: ";
    for (float value : camTemp) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    std::cout << "CpuTemp: " << event.GetCpuTemp() << std::endl;
    std::cout << "CoboTemp: " << event.GetCoboTemp() << std::endl;
    std::cout << "RadTemp: " << event.GetRadTemp() << std::endl;
    std::cout << "lvpsVol: " << event.GetlvpsVol() << std::endl;
    std::cout << "pumpVol: " << event.GetpumpVol() << std::endl;
    std::cout << "lvpsCur: " << event.GetlvpsCur() << std::endl;
    std::cout << "pumpCur: " << event.GetpumpCur() << std::endl;
    const std::vector<int>& siabMPWR = event.GetsiabMPWR();
    std::cout << "siabMPWR: ";
    for (int value : siabMPWR) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    const std::vector<int>& hvSW = event.GethvSW();
    std::cout << "hvSW: ";
    for (int value : hvSW) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}
