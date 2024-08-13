#include "ICalibration.h"


ICalibration::ICalibration(const std::string& filename)
{
    //relativeGain = 0.0;
    //std::cout << filename << endl;
    std::ifstream file(filename);
    std::string line;
    
    // Skip the header line
    std::getline(file, line);

    // Read data lines and populate the vectors
    while (std::getline(file, line)) {
        std::istringstream iss(line);
       //std::cout << line << endl;
        double pixelID, amplitude, errAmplitude, charge, errCharge;
        iss >> pixelID >> amplitude >> errAmplitude >> charge >> errCharge;
        amplitudes.push_back(amplitude);
        errAmplitudes.push_back(errAmplitude);
        charges.push_back(charge);
        errCharges.push_back(errCharge);
    }
}

double ICalibration::GetAmplitude(int pixelID) const {
    return amplitudes[pixelID];
}

double ICalibration::GetErrAmplitude(int pixelID) const {
    return errAmplitudes[pixelID];
}

double ICalibration::GetCharge(int pixelID) const {
    return charges[pixelID];
}

double ICalibration::GetErrCharge(int pixelID) const {
    return errCharges[pixelID];
}

double ICalibration::GetRelativeGain(int pixelID, double sipmTemp, double OperatingVol) {
    OperatingVol = OperatingVol - VOffsetPx[pixelID]/1000.0; //mV to V
   double relativeGain = (1 / (1 - BrVperDeg*(sipmTemp - 25)/(OperatingVol - BrVol[pixelID/16])));
    return relativeGain;
}

double ICalibration::GetRelativeOverVoltage(int pixelID, double operatingVol, double nominalVol){

    return ((operatingVol - VOffsetPx[pixelID]/1000.0 - BrVol[pixelID/16])/(nominalVol - VOffsetPx[pixelID]/1000.0 - BrVol[pixelID/16]));
}

double ICalibration::GetOverVoltage(int pixelID, double operatingVol){
    return (operatingVol - VOffsetPx[pixelID]/1000.0 - BrVol[pixelID/16])/BrVol[pixelID/16];
}

double ICalibration::GetVOffsetPx(int pixelID){
    return VOffsetPx[pixelID];
}

std::vector<double> ICalibration::GetVOffsetPx(){
    int size = sizeof(VOffsetPx)/sizeof(VOffsetPx[0]);
    std::vector<double> volOffsets(VOffsetPx, VOffsetPx + size);
    return volOffsets;
}

double ICalibration::GetMeasurementFromFile(int pixelID, double sipmTemp, double OperatingVol, std::string filename, std::string measurement, std::string biasVoltageRef){

    std::transform(biasVoltageRef.begin(), biasVoltageRef.end(), biasVoltageRef.begin(), ::toupper);

    if(biasVoltageRef.compare("ABS")==0){
        OperatingVol = OperatingVol - VOffsetPx[pixelID]/1000.0;
    }else if(biasVoltageRef.compare("REL")==0){
        OperatingVol = GetOverVoltage(pixelID, OperatingVol)/GetRelativeGain(pixelID,sipmTemp,OperatingVol);
    }

    int closestTempLine = GetClosestTemperatureLine(sipmTemp,filename,measurement);

    //std::cout<<closestTempLine<<std::endl;
    int iLine = 0;
    std::ifstream file(filename);
    std::string line;

    bool isMeasurementFound = false;
    std::istringstream iss;
    double opVol, opVolErr, meas, measErr;
    double measuredValue;

    std::vector<double> vOpVol, vMeas;
    while(std::getline(file, line) && !isMeasurementFound){
        if(iLine == closestTempLine){
            std::getline(file, line);
            std::getline(file, line);
            std::getline(file, line);
            do{ 
                //std::cout<<line<<std::endl;
                iss = std::istringstream(line);
                iss >> opVol >> opVolErr >> meas >> measErr;
                vOpVol.push_back(opVol);
                vMeas.push_back(meas);
                std::getline(file, line);
            }while(!line.empty());
            isMeasurementFound = true;
        }
        iLine++;
    }

    int i = 0;

    while(OperatingVol > vOpVol[i] && i < vOpVol.size()){
        i++;
    }

    measuredValue = IUtilities::Interpolate(vOpVol[i-1],vMeas[i-1],vOpVol[i],vMeas[i],OperatingVol);
    file.close();
    return measuredValue;
}

int ICalibration::GetMeasurementLine(std::string filename, std::string measurement){
    std::transform(measurement.begin(), measurement.end(), measurement.begin(), ::toupper);

    std::ifstream file(filename);
    std::string line;
    int iLine = 0;
    int ctr = 0;
    do{
        std::getline(file,line);
       //std::cout<<line<<std::endl;
        iLine = ++ctr;
        std::transform(line.begin(),line.end(),line.begin(), ::toupper);
    }while(line.compare(measurement)!=0);
    
    file.close();

    return iLine;
}

int ICalibration::GetClosestTemperatureLine(double sipmTemp, std::string filename, std::string measurement){

    int iMeasurementLine = GetMeasurementLine(filename, measurement);
    //std::cout<<iMeasurementLine<<std::endl;
    std::ifstream file(filename);
    std::string line;
    double difference = 1000;
    int closestTemp = 0;
    bool isMeasurementFound = false;
    bool isMeasurementFinished = false;
    int iLine = 0;

    while(std::getline(file,line) && !isMeasurementFound){
        if(iLine == iMeasurementLine){
            isMeasurementFound = true;
            while(std::getline(file,line) && !isMeasurementFinished){
                iLine++;
                if(line == "# Temperature"){
                    std::getline(file, line);
                    iLine++;
                    if(std::abs(sipmTemp - std::stod(line)) < difference){
                        difference = std::abs(sipmTemp - std::stod(line));
                        closestTemp = iLine;
                    }
                }else if(line.compare("***") == 0){
                    isMeasurementFinished = true;
                }
            }
        }
        iLine++;
    }

    file.close();

    return closestTemp;
}


ICalibration::~ICalibration() {}