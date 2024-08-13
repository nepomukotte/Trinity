#include "IEvent.h"

IEvent::IEvent(){
	//Iazimuth = 0;
	//Ialtitude = 0;
	IcorrectionFactor = 0;

/*    timestampIndex = -1;
    timestamp2Index = -1;

    closestTimestamp1 = 0;
    closestTimestamp2 = 0;*/

    latitude = 0.0f;
    longitude = 0.0f;
    altitude = 0;
    sunAzimuth = 0.0f;
    sunElevation = 0.0f;
    moonAzimuth = 0.0f;
    moonElevation = 0.0f;
    horizon = 0.0f;
    azimuth = 0.0f;

    TrigEvent = 0;
    TempFlag = 0;
    emon1 = 0;
    emon2 = 0;
    hv.resize(8, 0.0f);
    hvc.resize(8, 0.0f);
    ucTemp.resize(32, 0.0f);
    CpuTemp = 0;
    CoboTemp = 0;
    RadTemp = 0;
    lvpsVol = 0.0f;
    pumpVol = 0.0f;
    lvpsCur = 0;
    pumpCur = 0;
    siabMPWR.resize(32, 0);
    hvSW.resize(32, 0);
    sipmTemp.resize(32, 0.0f);
}

// void IEvent::LoadFromFile(const std::string& filename1, const std::string& filename2) {
//     std::ifstream file1(filename1);
//     if (!file1.is_open()) {
//         throw std::runtime_error("Failed to open file: " + filename1);
//     }

//     std::string line;
//     if (std::getline(file1, line)) {
//         std::stringstream ss(line);
//         std::string header;
//         while (ss >> header) {
//             headers.push_back(header);
//         }
//     }

//     if (headers.empty() || headers.size() != 9) {
//         throw std::runtime_error("Invalid file format: " + filename1);
//     }

//     timestampIndex = 0;
//     for (int i = 0; i < headers.size(); ++i) {
//         if (headers[i] == "timestamp") {
//             timestampIndex = i;
//             break;
//         }
//     }

//     while (std::getline(file1, line)) {
//         std::stringstream ss(line);
//         std::string cell;
//         std::vector<std::string> row;

//         while (ss >> cell) {
//             row.push_back(cell);
//         }

//         if (row.size() != headers.size()) {
//             throw std::runtime_error("Invalid file format: " + filename1);
//         }

//         data.push_back(row);
//     }

//     file1.close();

//     std::ifstream file2(filename2);
//     if (!file2.is_open()) {
//         throw std::runtime_error("Failed to open file: " + filename2);
//     }

//     if (std::getline(file2, line)) {
//         std::stringstream ss(line);
//         std::string header;
//         while (ss >> header) {
//             headers2.push_back(header);
//         }
//     }

//     if (headers2.empty() || headers2.size() != 124) {
//         throw std::runtime_error("Invalid file format: " + filename2);
//     }

//     timestamp2Index = 0;
//     for (int i = 0; i < headers2.size(); ++i) {
//         if (headers2[i] == "timestamp2") {
//             timestamp2Index = i;
//             break;
//         }
//     }

//     while (std::getline(file2, line)) {
//         std::stringstream ss(line);
//         std::string cell;
//         std::vector<std::string> row;

//         while (ss >> cell) {
//             row.push_back(cell);
//         }

//         if (row.size() != headers2.size()) {
//             throw std::runtime_error("Invalid file format: " + filename2);
//         }

//         data2.push_back(row);
//     }

//     file2.close();
// }

void IEvent::SetParametersFromTimestamp(int timestamp,std::vector<std::vector<std::string>> data, std::vector<std::vector<std::string>> data2) {
    // Find closest timestamp in the first file
    int closestTimestamp1 = findClosestTimestamp(data, timestamp);
   /* if (closestTimestamp1 != timestamp) {
        std::cout << "Given timestamp not found in the first file. Closest timestamp is: " << closestTimestamp1 << std::endl;
    }*/

    // Find closest timestamp in the second file
    int closestTimestamp2 = findClosestTimestamp(data2, timestamp);
    /*if (closestTimestamp2 != timestamp) {
        std::cout << "Given timestamp not found in the second file. Closest timestamp is: " << closestTimestamp2 << std::endl;
    }*/

    // Set parameters from the first file
    for (int i = 0; i < data.size(); ++i) {
        if (std::stoi(data[i][0]) == closestTimestamp1) {
            latitude = std::stof(data[i][1]);
            longitude = std::stof(data[i][2]);
            altitude = std::stoi(data[i][3]);
            sunAzimuth = std::stof(data[i][4]);
            sunElevation = std::stof(data[i][5]);
            moonAzimuth = std::stof(data[i][6]);
            moonElevation = std::stof(data[i][7]);
            horizon = std::stof(data[i][8]);
            azimuth = std::stof(data[i][9]);
            break;
        }
    }

    // Set parameters from the second file
    for (int i = 0; i < data2.size(); ++i) {
        if (std::stoi(data2[i][0]) == closestTimestamp2) {
            TrigEvent = std::stoi(data2[i][1]);
            TempFlag = std::stoi(data2[i][2]);
            emon1 = std::stoi(data2[i][3]);
            emon2 = std::stoi(data2[i][4]);
            for (int j = 0; j < 8; ++j) {
                hv[j] = std::stof(data2[i][j + 5]);
            }
            for (int j = 0; j < 8; ++j) {
                hvc[j] = std::stof(data2[i][j + 13]);
            }
            for (int j = 0; j < 32; ++j) {
                ucTemp[j] = std::stof(data2[i][j + 21]);
            }
            CpuTemp = std::stoi(data2[i][53]);
            CoboTemp = std::stoi(data2[i][54]);
            RadTemp = std::stoi(data2[i][55]);
            lvpsVol = std::stof(data2[i][56]);
            pumpVol = std::stof(data2[i][58]);
            lvpsCur = std::stoi(data2[i][57]);
            pumpCur = std::stoi(data2[i][59]);
            for (int j = 0; j < 32; ++j) {
                siabMPWR[j] = std::stoi(data2[i][j + 60]);
            }
            for (int j = 0; j < 32; ++j) {
                hvSW[j] = std::stoi(data2[i][j + 91]);
            }
            for (int j = 0; j < 32; ++j) {
                sipmTemp[j] = std::stof(data2[i][j + 122]);
            }
            break;
        }
    }
}

void IEvent::SetTelescopePointing (int timestamp,std::vector<std::vector<std::string>> data){

    int nearTimeStamp = findClosestTimestamp(data, timestamp);
    bool isFound = false;
    int i = 0;

   while(!isFound && i<data.size()){
        if (std::stoi(data[i][0]) == nearTimeStamp) {
            tiltAngle = std::stof(data[i][3]);
            isFound = true;
        }
        i++;
    }

}

void IEvent::SetTelescopePointingRaw(int timestamp,std::vector<std::vector<std::string>> data){

    int nearTimeStamp = findClosestTimestamp(data, timestamp);
    bool isFound = false;
    int i = 0;

   while(!isFound && i<data.size()){
        if (std::stoi(data[i][0]) == nearTimeStamp) {
            tiltAngleRaw = std::stof(data[i][1]);
            isFound = true;
        }
        i++;
    }

}

int IEvent::GetclosestTimestamp1() const {
    return 0;
}

int IEvent::GetclosestTimestamp2() const {
    return 0;
}

float IEvent::GetLatitude() const {
    return latitude;
}

float IEvent::GetLongitude() const {
    return longitude;
}

int IEvent::GetAltitude() const {
    return altitude;
}

float IEvent::GetSunAzimuth() const {
    return sunAzimuth;
}

float IEvent::GetSunElevation() const {
    return sunElevation;
}

float IEvent::GetMoonAzimuth() const {
    return moonAzimuth;
}

float IEvent::GetMoonElevation() const {
    return moonElevation;
}

float IEvent::GetHorizon() const {
    return horizon;
}

float IEvent::GetAzimuth() const {
    return azimuth;
}

int IEvent::GetTrigEvent() const {
    return TrigEvent;
}

int IEvent::GetTempFlag() const {
    return TempFlag;
}

int IEvent::GetEmon1() const {
    return emon1;
}

int IEvent::GetEmon2() const {
    return emon2;
}

const std::vector<float>& IEvent::Gethv() const {
    return hv;
}

const std::vector<float>& IEvent::Gethvc() const {
    return hvc;
}

const std::vector<float>& IEvent::GetUCTemp() const {
    return ucTemp;
}

int IEvent::GetCpuTemp() const {
    return CpuTemp;
}

int IEvent::GetCoboTemp() const {
    return CoboTemp;
}

int IEvent::GetRadTemp() const {
    return RadTemp;
}

float IEvent::GetlvpsVol() const {
    return lvpsVol;
}

float IEvent::GetpumpVol() const {
    return pumpVol;
}

int IEvent::GetlvpsCur() const {
    return lvpsCur;
}

int IEvent::GetpumpCur() const {
    return pumpCur;
}

const std::vector<int>& IEvent::GetsiabMPWR() const {
    return siabMPWR;
}

const std::vector<int>& IEvent::GethvSW() const {
    return hvSW;
}

const std::vector<float>& IEvent::GetSiPMTemp() const {
    return sipmTemp;
}

int IEvent::findClosestTimestamp(const std::vector<std::vector<std::string>>& data, int timestamp) const {
    int closestTimestamp = std::stoi(data[0][0]);
    int minDifference = std::abs(timestamp - closestTimestamp);

    for (int i = 1; i < data.size(); ++i) {
        int currentTimestamp = std::stoi(data[i][0]);
        int difference = std::abs(timestamp - currentTimestamp);

        if (difference < minDifference) {
            minDifference = difference;
            closestTimestamp = currentTimestamp;
        }
    }

    return closestTimestamp;
}

int IEvent::FindClosestTimestampIndex(std::vector<unsigned long long> timeArray, unsigned long long eventTime){
    unsigned long long closestTimestamp = timeArray[0];
    unsigned long long minDifference = std::abs((long long int)(eventTime - closestTimestamp));

    unsigned long long currentTimestamp;
    unsigned long long difference;

    int index = 0;
    for (int i = 1; i < timeArray.size(); ++i) {
        currentTimestamp = timeArray[i];
        
        difference = std::abs((long long int)(eventTime - currentTimestamp));

        if (difference < minDifference) {
            minDifference = difference;
            closestTimestamp = currentTimestamp;
            index = i;
        }
    }

    return index;
}

IEvent::~IEvent(){

}

float IEvent::GetTiltAngle() const{
    return tiltAngle;
}

float IEvent::GetTiltAngleRaw() const{
    return tiltAngleRaw;
}