#include <algorithm>
#include "commands.h"
#include "fcutils/utility/utility.hpp"
#include "fcutils/string/string.hpp"
#include "smbus.c"

//int file = 0;
int I2C_ADDR = 0;
int First_SIAB_ID = 0x10;
int SIAB_ID_List[32] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,32,31};

int CMD_Format_1[] = {0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x10,0x11,0x12,
                      0x13,0x15,0x16,0x17,0x18,0x19,0x1C,0x1E,0x1F,0x20,0x22,0x23,0x24};
int CMD_Format_2[] = {0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0x14,0x1A,0x21};
int CMD_Format_3[] = {0x1B,0x1D};

size_t format_size_1 = sizeof(CMD_Format_1)/sizeof(int);
size_t format_size_2 = sizeof(CMD_Format_2)/sizeof(int);
size_t format_size_3 = sizeof(CMD_Format_3)/sizeof(int);

using namespace std;

void Reset_Atmega328P(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "Reset command is sent successfully." << std::endl;
    }else{
        std::cout << "Sending Reset command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Turn_on_LED(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "LED Turn On command is sent successfully." << std::endl;
    }else{
        std::cout << "Sending LED Turn On command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Turn_off_LED(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "LED Turn Off command is sent successfully." << std::endl;
    }else{
        std::cout << "Sending LED Turn Off command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Turn_on_3V3(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "Music 3.3V On command is sent successfully." << std::endl;
    }else{
        std::cout << "Sending Music 3.3v On command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Turn_off_3V3(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "Music 3.3V Off command is sent successfully." << std::endl;
    }else{
        std::cout << "Sending Music 3.3v Off command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Turn_on_5V(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "Music 5V On command is sent successfully." << std::endl;
    }else{
        std::cout << "Sending Music 5v On command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Turn_off_5V(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "Music 5V Off command is sent successfully." << std::endl;
    }else{
        std::cout << "Sending Music 5v Off command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Turn_on_3V3_5V(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "Music 3.3V and 5V On command is sent successfully." << std::endl;
    }else{
        std::cout << "Sending Music 3.3V and 5V On command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Turn_off_3V3_5V(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "Music 3.3V and 5V Off command is sent successfully." << std::endl;
    }else{
        std::cout << "Sending Music 3.3V and 5V Off command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Read_SiPM_Temp(int file, __u32 command, std::string &response)
{
    int status = 0;
    int temperature = 0;
    __u8 cmd_code = command >> 24;
    __s32 result = i2c_smbus_read_word_data(file, cmd_code);
    if(result < 0){
        status = result;
        temperature = -1;
        std::cout << "Read SiPM Temperature command failed." << std::endl;
    }else{
        status = 0;
        temperature = result;
        std::cout << "SiPM Temp: " << temperature << std::endl;
    }

    FILE *filePointer;
    char filename[100] = "sipm_tmp.txt";
    filePointer = fopen(filename, "a");
    fwrite(&temperature, sizeof(int), 1, filePointer);
    fclose(filePointer);

    auto temp = fcutils::string::int_to_hex(status) + fcutils::string::int_to_hex(temperature); 
    //auto temp = fcutils::string::int_to_hex(temperature);
    response += temp;
}

void Read_UC_Temp(int file, __u32 command, std::string &response)
{
    int status = 0;
    int temperature = 0;
    __u8 cmd_code = command >> 24;
    __s32 result = i2c_smbus_read_word_data(file, cmd_code);
    if(result < 0){
        status = result;
        temperature = -1;
        std::cout << "Read UC Temperature command failed." << std::endl;
    }else{
        status = 0;
        temperature = result;
        std::cout << "UC Temp: " << temperature << std::endl;
    }

    FILE *filePointer;
    char filename[100] = "uc_tmp.txt";
    filePointer = fopen(filename, "a");
    fwrite(&temperature, sizeof(int), 1, filePointer);
    fclose(filePointer);

    auto temp = fcutils::string::int_to_hex(status) + fcutils::string::int_to_hex(temperature);
    //auto temp = fcutils::string::int_to_hex(temperature);
    response += temp;
}

__u8 Read_EPRM_Byte(int file, __u32 command, std::string &response)
{
    int status = 0;
    int value = 0;
    __u8 cmd_code = ((command >> 24) & 0xFF);
    __u16 cmd_data = ((command >> 8) & 0xFFFF);
    __s32 result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    if(result < 0){
        status = result;
        value = -1;
        std::cout << "Read EPRM Byte command failed." << std::endl;
    }else{
        status = 0;
        value = result;
        std::cout << "Read one byte from EEPROM at address " << cmd_data << ": " << value << std::endl;    
    }
    auto temp = fcutils::string::int_to_hex(status) + fcutils::string::int_to_hex(value);
    response += temp;

    return (value & 0xFF);
}

__u16 Read_EPRM_Word(int file, __u32 command, std::string &response)
{
    int status = 0;
    int value = 0;
    __u8 cmd_code = ((command >> 24) & 0xFF);
    __u16 cmd_data = ((command >> 8) & 0xFFFF);
    __s32 result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    if(result < 0){
        status = result;
        value = -1;
        std::cout << "Read EPRM Word command failed." << std::endl;
    }else{
        status = 0;
        value = result;
        std::cout << "Read one word from EEPROM at address " << cmd_data << ": " << value << std::endl;    
    }
    auto temp = fcutils::string::int_to_hex(status) + fcutils::string::int_to_hex(value);
    response += temp;
    
    return (result & 0xFFFF);
}

// For writing one byte, 2nd and 3rd byte of command are the address and 4th byte is the byte value.
void Write_EPRM_Byte(int file, __u32 command, std::string &response)
{
    __u8 cmd_data[3];
    __u8 cmd_code = ((command >> 24) & 0xFF);
    cmd_data[0] = ((command >> 8) & 0xFF);      // this is the low byte of address we are writing to 
    cmd_data[1] = ((command >> 16) & 0xFF);     // this is the high byte of address we are writing to    
    cmd_data[2] = (command & 0xFF);             // this is the byte value we are writing   
    int status = i2c_smbus_write_block_data(file, cmd_code, 3, cmd_data);
    if(status == 0){
        std::cout << "One byte is written to EEPROM." << std::endl;
    }else{
        std::cout << "Writing one byte to EEPROM command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

// For writing one word, 2nd byte of command is the address and 3rd and 4th are word value.
// Due to limitation on command length, we use only 256 bytes of EEPROM for now.
void Write_EPRM_Word(int file, __u32 command, std::string &response)
{
    __u8 cmd_data[3];
    __u8 cmd_code = ((command >> 24) & 0XFF);
    cmd_data[0] = ((command >> 16) & 0xFF);     // this is the address we are writing to. Note that the address should be only one byte here.
    cmd_data[1] = (command & 0xFF);             // this is the low byte of value we are writing
    cmd_data[2] = ((command >> 8) & 0xFF);      // this is the high byte of value we are writing
    int status = i2c_smbus_write_block_data(file, cmd_code, 3, cmd_data);
    if(status == 0){
        std::cout << "One word is written to EEPROM." << std::endl;
    }else{
        std::cout << "Writing one word to EEPROM command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Init_Music(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = ((command >> 24) & 0XFF);
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "Music SPI is initialized successfully." << std::endl;
    }else{
        std::cout << "Music SPI initialization command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Load_Music_Config_EPRM(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = ((command >> 24) & 0XFF);
    __u8 cmd_data = ((command >> 16) & 0xFF);
    int status = i2c_smbus_write_byte_data(file, cmd_code, cmd_data);
    if (status == 0 && cmd_data == 0x01){
        std::cout << "Bottom Music Configuration is loaded from EEPROM." << std::endl;
    }else if (status == 0 && cmd_data == 0x02){
        std::cout << "Top Music Configuration is loaded from EEPROM." << std::endl;
    }else{
        std::cout << "An error occured. Loading configuration from EEPROM failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Print_Loaded_Music_Config(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = ((command >> 24) & 0XFF);
    __u8 cmd_data = ((command >> 16) & 0xFF);
    int status = i2c_smbus_write_byte_data(file, cmd_code, cmd_data);
    if (status == 0 && cmd_data == 0x01){
        std::cout << "Bottom Music Configuration is sent to slave output." << std::endl;
    }else if (status == 0 && cmd_data == 0x02){
        std::cout << "Top Music Configuration is sent to slave output." << std::endl;
    }else{
        std::cout << "An error occured. Loading configuration from EEPROM failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Write_Loaded_Config_to_Music(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = ((command >> 24) & 0XFF);
    __u8 cmd_data = ((command >> 16) & 0xFF);
    int status = i2c_smbus_write_byte_data(file, cmd_code, cmd_data);
    if (status == 0 && cmd_data == 0x01){
        std::cout << "Configuration loaded from EEPROM is written to bottom Music." << std::endl;
    }else if (status == 0 && cmd_data == 0x02){
        std::cout << "Configuration loaded from EEPROM is written to top Music." << std::endl;
    }else{
        std::cout << "An error occured. Writing configuration to Music from EEPROM failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

__u16 Read_from_Music_Register(int file, __u32 command, std::string &response)
{
    int status;
    int value;
    __u8 cmd_code = ((command >> 24) & 0XFF);
    __u16 cmd_data = ((command >> 8) & 0xFFFF);
    __s32 result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    if(result < 0){
        status = result;
        value = -1;
        std::cout << "Read from Music Register command failed." << std::endl;
    }else if((result >= 0) && (((cmd_data >> 8) & 0xFF) == 0x01)){
        status = 0;
        value = result;
        std::cout << "Reading from Bottom Music Register: " << value << std::endl;
    }else if((result >= 0) && (((cmd_data >> 8) & 0xFF) == 0x02)){
        status = 0;
        value = result;
        std::cout << "Reading from Top Music Register: " << value << std::endl;
    }else{
        status = -1;
        std::cout << "An error occured. Unknown state!" << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status) + fcutils::string::int_to_hex(value);
    response += temp;

    return value;
}

bool Write_to_Music_Register(int file, __u32 command, std::string &response)
{
    __u8 cmd_data[3];
    __u8 cmd_code = ((command >> 24) & 0xFF);
    cmd_data[0] = ((command >> 16) & 0xFF);
    cmd_data[1] = (command & 0xFF);             // This is the LSB of value we are writing to Music
    cmd_data[2] = ((command >> 8) & 0xFF);      // This is the MSB of value we are writing to Music
    int status = i2c_smbus_write_block_data(file, cmd_code, 3, cmd_data);
    if ((status == 0) && ((cmd_data[0] >> 6) == 0x01)){
        std::cout << "Writing to bottom Music register was successful." << std::endl;
    }else if ((status == 0) && ((cmd_data[0] >> 6) == 0x02)){
        std::cout << "Writing to top Music register was successful." << std::endl;
    }else{
        std::cout << "An error occured. Writing to Music register failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;

    if (status == 0){
        return true;
    }else{
        return false;
    }
}

void Set_ADC_Rate(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = ((command >> 24) & 0XFF);
    __u8 cmd_data = ((command >> 16) & 0xFF);
    int status = i2c_smbus_write_byte_data(file, cmd_code, cmd_data);
    if(status == 0){
        std::cout << "ADC rate is set successfully." << std::endl;
    }else{
        std::cout << "Sending ADC rate set command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Initialize_ADC(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "ADC is initialized successfully." << std::endl;
    }else{
        std::cout << "Sending ADC initialization command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Enable_ADC_Run(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "ADC Interrupt Service Routine is activated." << std::endl;
    }else{
        std::cout << "Sending ADC Interrupt Service Routine command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Disable_ADC_Run(int file, __u32 command, std::string &response)
{
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "ADC Interrupt Service Routine is disabled." << std::endl;
    }else{
        std::cout << "Disabling ADC Interrupt Service Routine command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Read_Single_Current(int file, __u32 command, std::string &response)
{
    int value;
    int status;
    __u8 cmd_code = ((command >> 24) & 0xFF);
    __u8 cmd_data = ((command >> 16) & 0xFF);
    __s32 result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    if(result < 0){
        status = result;
        value = -1;
        std::cout << "Reading single current command failed." << std::endl;
    }else{
        status = 0;
        value = ((result << 8) | 0x00);
        std::cout << "Reading single current was successful." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status) + fcutils::string::int_to_hex(value);
    response += temp;

    double converted_value = value * 0.000149011612;
    std::cout << "ADC Sample from channel " << cmd_data << " is received: " << converted_value << std::endl;
}

void Read_All_Current(int file, __u32 command, std::string &response)
{
    __u32 cmd_tmp;
    __u8 cmd_code = ((command >> 24) & 0xFF);
    for (int i=0; i < 16; i++){
        cmd_tmp = ((CMD_READ_SINGLE_CURRENT << 24) | (i << 16) | 0x0000);
        Read_Single_Current(file, cmd_tmp, response);
        cmd_tmp=0;
    }
}

void Record_ADC_Offset(int file, __u32 command, std::string &response)
{
    __s32 result;
    __u8 cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "ADC offset is recorded successfully." << std::endl;
    }else{
        std::cout << "ADC offset record command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Retrieve_ADC_Offset(int file, __u32 command, std::string &response)
{
    int offset;
    int status;
    __u8 cmd_data;
    __s32 result;
    __u8 cmd_code = ((command >> 24) & 0xFF);
    for (__u8 i=0; i < 16; i++){
        cmd_data = i;
        result = i2c_smbus_process_call(file, cmd_code, cmd_data);
        if(result < 0){
            status = -1;
            offset = -1;
            std::cout << "Retrieve ADC offset command from channel#" << unsigned(cmd_data) << " failed." << std::endl;
        }else{
            status = 0;
            offset = ((result << 8) | 0x00);
            std::cout << "ADC Offset from channel#" << unsigned(cmd_data) << " is received: " << double(offset) * 0.000149011612 << std::endl;
        }
        auto temp = fcutils::string::int_to_hex(status) + fcutils::string::int_to_hex(offset);
        response += temp;
    }
}

void Turn_on_HV(int file, __u32 command, std::string &response)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "HV Switch on is sent successfully." << std::endl;
    }else{
        std::cout << "HV Switch on command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Turn_off_HV(int file, __u32 command, std::string &response)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "HV Switch off is sent successfully." << std::endl;
    }else{
        std::cout << "HV Switch off command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Enable_HV_CTR(int file, __u32 command, std::string &response)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    int status = i2c_smbus_write_byte(file, cmd_code);
    if(status == 0){
        std::cout << "Enable HV Control is sent successfully." << std::endl;
    }else{
        std::cout << "Enable HV Control command failed." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

__u8 Check_HV_Status(int file, __u32 command, std::string &response)
{
    int status;
    int value;
    __u8 cmd_code;
    __u8 HV_Status;
    __s32 result;
    cmd_code = ((command >> 24) & 0xFF);
    result = i2c_smbus_process_call(file, cmd_code, 0x00);
    if(result < 0){
        status = result;
        value = 0;
        std::cout << "Check HV status command failed." << std::endl;
    }else if(result == 1){
        status = 0;
        value = result;
        std::cout << "HV is off." << std::endl;
    }else if(result == 2){
        status = 0;
        value = result;
        std::cout << "HV is on." << std::endl;
    }else{
        status = -1;
        value = -1;
        std::cout << "Unknown State." << std::endl;
    }

    FILE *filePointer;
    char filename[100] = "hv_state_tmp.txt";
    filePointer = fopen(filename, "a");
    fwrite(&value, sizeof(int), 1, filePointer);
    fclose(filePointer);

    auto temp = fcutils::string::int_to_hex(status) + fcutils::string::int_to_hex(value);
    response += temp;

    return result;
}

void Set_Max_Current(int file, __u32 command, std::string &response)
{
    int status;
    __s32 result;
    __u8 cmd_code;
    __u16 cmd_data;
    cmd_code = ((command >> 24) & 0xFF);
    cmd_data = ((command >> 8) & 0xFFFF);
    result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    if(result < 0){
        status = result;
        std::cout << "Set max current command failed." << std::endl;
    }else{
        status = 0;
        std::cout << "Set max current command sent successfully." << std::endl;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

bool Read_Music_Config(std::string ConfigCSV, int *address, int *value)
{
    std::ifstream MusicConfig(ConfigCSV);
    if(!MusicConfig.is_open()){
        std::cout << "Can not open Music Configuration file." << std::endl;
        return false;
    }

    std::string line, address_str, value_str;

    // Reading the first 8 register addresses and values
    getline(MusicConfig,line);
    for(int i=0; i<8; i++)
    {
        getline(MusicConfig,line);
        address_str = line.substr(0,line.find(","));
        line = line.substr(line.find(",")+1);
        value_str = line.substr(0,line.find(","));

        address[i] = stoi(address_str);
        value[i] = stoi(value_str);
    }

    // Reading the second 8 register addresses and values
    getline(MusicConfig,line);
    for(int i=8; i<16; i++)
    {
        getline(MusicConfig,line);
        address_str = line.substr(0,line.find(","));
        line = line.substr(line.find(",")+1);
        value_str = line.substr(0,line.find(","));

        address[i] = stoi(address_str);
        value[i] = stoi(value_str);
    }

    // Reading the last 7 register addresses and values
    for(int i=16; i<23; i++)
    {
        getline(MusicConfig,line);
        getline(MusicConfig,line);
        address_str = line.substr(0,line.find(","));
        line = line.substr(line.find(",")+1);
        value_str = line.substr(0,line.find(","));

        address[i] = stoi(address_str);
        value[i] = stoi(value_str);
    }

    MusicConfig.close();
    return true;
}

void Update_Music_Config(int file, __u32 command, std::string &response)
{
    int address[23];
    int value[23];
    int write_cmd = 21;
    int fail_ctr = 0;
    int status;
    __u32 cmd_rebuilt;
    std::string response_fake;    // this is being passed, so the Write_to_Music_Register function does not complain, but we are not using it.
    std::string filename= "/home/cherenkov/Programs/control_software/siab/Music_Register_Table.csv";
    if(!Read_Music_Config(filename, address, value)){
        status = -60;
    }else{
        __u8 cmd_code = ((command >> 24) & 0xFF);
        __u8 music_ID = ((command >> 16) & 0xFF);

        for (int i=0; i<23; i++)
        {
            int value_LSB = (value[i] & 0xFF);
            int value_MSB = ((value[i] >> 8) & 0xFF);
            cmd_rebuilt = (((__u32)write_cmd << 24) | ((__u32)music_ID << 22) | ((__u32)address[i] << 16) | ((__u32)value_MSB << 8) | (__u32)(value_LSB));
            if(!Write_to_Music_Register(file, cmd_rebuilt, response_fake))
                fail_ctr++;
            usleep(100000);
        }
        status = fail_ctr;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Set_Bias_Voltage(int i, int file, __u32 command, std::string &response)
{
    int status = 0;
    std::string filename = "/home/cherenkov/Programs/control_software/siab/DV_Offset.csv";
    std::ifstream BiasVoltage(filename);
    if(!BiasVoltage.is_open()){
        status = -60;           // this is just to make sure error code for file opening failure is different from other error codes.
        std::cout << "Can not open Bias Voltage Configuration file." << std::endl;
    }else{
        int fail_ctr = 0;
        int startline = i*16;
        int BV_value[16] = {0};
        std::string BV_str;

        std::string response_fake;    // this is being passed, so the Write_to_Music_Register function does not complain, but we are not using it.

        // Reading the 16 Bias voltage value for the SIAB pointed by i
        getline(BiasVoltage,BV_str);
        for(int j=0; j<startline; j++)
        {
            getline(BiasVoltage,BV_str);
        }
        for(int k=startline; k<startline+16; k++)
        {
            getline(BiasVoltage,BV_str);
            BV_value[k-startline] = stoi(BV_str);
        }
        BiasVoltage.close();

        std::cout << std::endl << "Setting bias voltage for SIAB# " << i << std::endl;

        __u32 cmd_rebuilt = 0;
        for (int k=0; k<8; k++)
        {
            cmd_rebuilt = (((__u32)0x15 << 24) | ((__u32)(0x01) << 22) | ((__u32)k << 16) | ((__u32)BV_value[k]));
            std::cout << std::hex << cmd_rebuilt << "\t";
            if(!Write_to_Music_Register(file, cmd_rebuilt, response_fake))
                fail_ctr++;
            usleep(100000);
            cmd_rebuilt = (((__u32)0x15 << 24) | ((__u32)(0x02) << 22) | ((__u32)k << 16) | ((__u32)BV_value[k+8]));
            std::cout << std::hex << cmd_rebuilt << std::endl;
            if(Write_to_Music_Register(file, cmd_rebuilt, response_fake))
                fail_ctr++;
            usleep(100000);
        }
        std::cout << "Finished Setting bias voltage for SIAB# " << i << std::endl;
        if(fail_ctr != 0)
            std::cout << "Bias voltage setup failed for " << fail_ctr << "channels on above SIAB." << std::endl;
        status = -fail_ctr;
    }
    auto temp = fcutils::string::int_to_hex(status);
    response += temp;
}

void Set_Trigger_Threshold(int file, __u32 command, std::string &response)
{
    __u32 CMD_Packet = 0;
    __u8 music_ID = ((command >> 16) & 0xFF);
    __u16 threshold_value = (command & 0xFFFF);
    threshold_value = ((threshold_value << 3) | 0xF007);
    std::string response_fake;    // this is being passed, so the Write_to_Music_Register function does not complain, but we are not using it.
    if(music_ID == 1){
        for (int j=8; j<16; j++){
            CMD_Packet = ((0x15 << 24) | (((0x01 << 6) | j) << 16) | threshold_value);
            Write_to_Music_Register(file, CMD_Packet, response_fake);
            usleep(10000);
        }
    }else if(music_ID == 2){
        for (int j=8; j<16; j++){
            CMD_Packet = ((0x15 << 24) | (((0x02 << 6) | j) << 16) | threshold_value);
            Write_to_Music_Register(file, CMD_Packet, response_fake);
            usleep(10000);
        }
    }else if(music_ID == 3){
        for (int j=8; j<16; j++){
            CMD_Packet = ((0x15 << 24) | (((0x02 << 6) | j) << 16) | threshold_value);
            Write_to_Music_Register(file, CMD_Packet, response_fake);
            usleep(10000);
            CMD_Packet = 0;
            CMD_Packet = ((0x15 << 24) | (((0x01 << 6) | j) << 16) | threshold_value);
            Write_to_Music_Register(file, CMD_Packet, response_fake);
            usleep(10000);
        }
    }else{
        std::cout << "An error occured. Could not understand the Music IDs." << std::endl;
    }
}

int Get_CMD_Format(int CMD_ID)
{
    int *end, *result;
    int CMD_Format = -1;
    bool found_flag = false;

    while(!found_flag){
        end = CMD_Format_1 + format_size_1;
        result = std::find(CMD_Format_1, end, CMD_ID);
        if(result != end){
            found_flag = true;
            CMD_Format = 1;
            break;
        }else{
            end = CMD_Format_2 + format_size_2;
            result = std::find(CMD_Format_2, end, CMD_ID);
            if(result != end){
                found_flag = true;
                CMD_Format = 2;
                break;
            }else{
                end = CMD_Format_3 + format_size_3;
                result = std::find(CMD_Format_3, end, CMD_ID);
                if(result != end){
                    found_flag = true;
                    CMD_Format = 3;
                    break;
                }else{
                    CMD_Format = 0;
                    break;
                }
            }
        }
    }
    return CMD_Format;
}

void Process_CMD(int file, std::vector<uint32_t>& vCMD, std::string &response)
{
    uint32_t SIAB_ID = ((vCMD.at(0) << 24) | (vCMD.at(1) << 16) | (vCMD.at(2) << 8) | (vCMD.at(3)));
    uint32_t CMD_ID = vCMD.at(4);
    uint32_t CMD_Packet = ((vCMD.at(4) << 24) | (vCMD.at(5) << 16) | (vCMD.at(6) << 8) | (vCMD.at(7)));

    int CMD_Format = Get_CMD_Format(CMD_ID);
/*
    auto response_id = fcutils::string::int_to_hex(CMD_Format);
    auto hex_id = fcutils::string::int_to_hex(SIAB_ID);
    response = response_id + hex_id + response;
*/
    for (int i=0; i<32; i++)
    {
        int SIAB_bit = 0;
        SIAB_bit = ((SIAB_ID >> i) & 0x01);
        if (SIAB_bit)
        {
            usleep(30000);
            std::string temp_response = "0x20202020";
            I2C_ADDR = First_SIAB_ID + SIAB_ID_List[i];
            if (ioctl(file, I2C_SLAVE, I2C_ADDR) < 0)
            {
                /*
                if(CMD_Format == 1){
                    temp_response = "0x" + std::string(2,'9').append(temp_response);
                }else if(CMD_Format == 2){
                    temp_response = "0x" + std::string(8,'6').append(temp_response);
                    temp_response = "0x" + std::string(2,'9').append(temp_response);
                } else if(CMD_Format == 3){
                    for(int j=0; j<16; j++){
                    temp_response = "0x" + std::string(8,'6').append(temp_response);
                    temp_response = "0x" + std::string(2,'9').append(temp_response);
                    }
                }else{
                    temp_response = "";
                }
                */
                response = response + "0x00000000";
                std::cout << "Error in setting up SIAB Number: " << i+1 << std::endl;
            }
            else{
                switch (CMD_ID)
                {
                    case CMD_RESET:
                        Reset_Atmega328P(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_ON_LED:
                        Turn_on_LED(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_OFF_LED:
                        Turn_off_LED(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_ON_3V3:
                        Turn_on_3V3(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_OFF_3V3:
                        Turn_off_3V3(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_ON_5V:
                        Turn_on_5V(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_OFF_5V:
                        Turn_off_5V(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_ON_3V3_5V:
                        Turn_on_3V3_5V(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_OFF_3V3_5V:
                        Turn_off_3V3_5V(file, CMD_Packet, temp_response);
                        break;
                    case CMD_READ_SIPM_TEMP:
                        Read_SiPM_Temp(file, CMD_Packet, response);
                        break;
                    case CMD_READ_UC_TEMP:
                        Read_UC_Temp(file, CMD_Packet, response);
                        break;
                    case CMD_READ_EPRM_BYTE:
                        Read_EPRM_Byte(file, CMD_Packet, temp_response);
                        break;
                    case CMD_READ_EPRM_WORD:
                        Read_EPRM_Word(file, CMD_Packet, temp_response);
                        break;
                    case CMD_WRITE_EPRM_BYTE:
                        Write_EPRM_Byte(file, CMD_Packet, temp_response);
                        break;
                    case CMD_WRITE_EPRM_WORD:
                        Write_EPRM_Word(file, CMD_Packet, temp_response);
                        break;
                    case CMD_INIT_MUSIC:
                        Init_Music(file, CMD_Packet, temp_response);
                        break;
                    case CMD_LOAD_MUSIC_CONFIG_EPRM:
                        Load_Music_Config_EPRM(file, CMD_Packet, temp_response);
                        break;
                    case CMD_PRINT_LOADED_MUSIC_CONFIG:
                        Print_Loaded_Music_Config(file, CMD_Packet, temp_response);
                        break;
                    case CMD_WRITE_LOADED_CONFIG_TO_MUSIC:
                        Write_Loaded_Config_to_Music(file, CMD_Packet, temp_response);
                        break;
                    case CMD_READ_FROM_MUSIC_REGISTER:
                        Read_from_Music_Register(file, CMD_Packet, temp_response);
                        break;
                    case CMD_WRITE_TO_MUSIC_REGISTER:
                        Write_to_Music_Register(file, CMD_Packet, temp_response);
                        break;
                    case CMD_SET_ADC_RATE:
                        Set_ADC_Rate(file, CMD_Packet, temp_response);
                        break;
                    case CMD_INITIALIZE_ADC:
                        Initialize_ADC(file, CMD_Packet, temp_response);
                        break;
                    case CMD_ENABLE_ADC_RUN:
                        Enable_ADC_Run(file, CMD_Packet, temp_response);
                        break;
                    case CMD_DISABLE_ADC_RUN:
                        Disable_ADC_Run(file, CMD_Packet, temp_response);
                        break;
                    case CMD_READ_SINGLE_CURRENT:
                        Read_Single_Current(file, CMD_Packet, temp_response);
                        break;
                    case CMD_READ_ALL_CURRENT:
                        Read_All_Current(file, CMD_Packet, temp_response);
                        break;
                    case CMD_RECORD_ADC_OFFSET:
                        Record_ADC_Offset(file, CMD_Packet, temp_response);
                        break;
                    case CMD_RETRIEVE_ADC_OFFSET:
                        Retrieve_ADC_Offset(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_ON_HV:
                        Turn_on_HV(file, CMD_Packet, temp_response);
                        break;
                    case CMD_TURN_OFF_HV:
                        Turn_off_HV(file, CMD_Packet, temp_response);
                        break;
                    case CMD_ENABLE_HV_CTR:
                        Enable_HV_CTR(file, CMD_Packet, temp_response);
                        break;
                    case CMD_CHECK_HV_STATUS:
                        Check_HV_Status(file, CMD_Packet, temp_response);
                        break;
                    case CMD_SET_MAX_CURRENT:
                        Set_Max_Current(file, CMD_Packet, temp_response);
                        break;
                    case CMD_UPDATE_MUSIC_CONFIG:
                        Update_Music_Config(file, CMD_Packet, temp_response);
                        break;
                    case CMD_SET_BIAS_VOLTAGE:
                        Set_Bias_Voltage(i, file, CMD_Packet, temp_response);
                        break;
                    case CMD_SET_TRG_THRESHOLD:
                        Set_Trigger_Threshold(file, CMD_Packet, temp_response);
                        break;
                    default:
                        break;
                }
            }
            //response = response + temp_response;
        }
    }
}
