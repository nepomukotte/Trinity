#include "commands.h"

void Reset_Atmega328P(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "Reset command is sent." << std::endl;
}

void Turn_on_LED(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "LED is on." << std::endl;    
}

void Turn_off_LED(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "LED is off." << std::endl;    
}

void Turn_on_3V3(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "Music 3.3V is on." << std::endl;    
}

void Turn_off_3V3(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "Music 3.3V is off." << std::endl;    
}

void Turn_on_5V(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "Music 5V is on." << std::endl;    
}

void Turn_off_5V(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "Music 5V is off." << std::endl;    
}

void Turn_on_3V3_5V(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "3.3V and 5V Switch is on." << std::endl;    
}

void Turn_off_3V3_5V(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "3.3V and 5V Switch is off." << std::endl;    
}

__u16 Read_SiPM_Temp(int file, __u32 command)
{
    __u8 cmd_code;
    __u16 tempeature;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_read_word_data(file, cmd_code);
    tempeature = result;
    std::cout << "SiPM Tempeature is: " << tempeature << std::endl;
    return tempeature;
}

__u16 Read_UC_Temp(int file, __u32 command)
{
    __u8 cmd_code;
    __u16 tempeature;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_read_word_data(file, cmd_code);
    tempeature = result;
    std::cout << "Microcontroller Tempeature is: " << tempeature << std::endl;
    return tempeature;
}

__u8 Read_EPRM_Byte(int file, __u32 command)
{
    __u8 cmd_code;
    __u16 cmd_data;
    __u8 ret_value;
    __s32 result;
    cmd_code = ((command >> 24) & 0xFF);
    cmd_data = ((command >> 8) & 0xFFFF);
    result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    ret_value = (result & 0xFF);
    std::cout << "Read one byte from EEPROM at address " << cmd_data << ": " << unsigned(ret_value) << std::endl;
    
    return result;
}

__u16 Read_EPRM_Word(int file, __u32 command)
{
    __u8 cmd_code;
    __u16 cmd_data;
    __u16 ret_value;
    __s32 result;
    cmd_code = ((command >> 24) & 0xFF);
    cmd_data = ((command >> 8) & 0xFFFF);
    result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    ret_value = (result & 0xFFFF);
    std::cout << "Read one word from EEPROM at address " << cmd_data << ": " << unsigned(ret_value) << std::endl;
    
    return result;
}

// For writing one byte, 2nd and 3rd byte of command are the address and 4th byte is the byte value.
void Write_EPRM_Byte(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 cmd_data[3];
    __s32 result;
    cmd_code = ((command >> 24) & 0xFF);
    cmd_data[0] = ((command >> 8) & 0xFF);      // this is the low byte of address we are writing to 
    cmd_data[1] = ((command >> 16) & 0xFF);     // this is the high byte of address we are writing to    
    cmd_data[2] = (command & 0xFF);             // this is the byte value we are writing   
    result = i2c_smbus_write_block_data(file, cmd_code, 3, cmd_data);
    if (result == 0)
    {
        std::cout << "One byte is written to EEPROM." << std::endl;
    }
}

// For writing one word, 2nd byte of command is the address and 3rd and 4th are word value.
// Due to limitation on command length, we use only 256 bytes of EEPROM for now.
void Write_EPRM_Word(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 cmd_data[3];
    __s32 result;
    cmd_code = ((command >> 24) & 0XFF);
    cmd_data[0] = ((command >> 16) & 0xFF);     // this is the address we are writing to. Note that the address should be only one byte here.
    cmd_data[1] = (command & 0xFF);             // this is the low byte of value we are writing
    cmd_data[2] = ((command >> 8) & 0xFF);      // this is the high byte of value we are writing
    std::cout << "cmd code: " << unsigned(cmd_code) << std::endl;
    std::cout << "Address: " << unsigned(cmd_data[0]) << std::endl;
    std::cout << "Word value: " << unsigned(command & 0xFFFF) << std::endl;
    result = i2c_smbus_write_block_data(file, cmd_code, 3, cmd_data);
    if (result == 0)
    {
        std::cout << "One word is written to EEPROM." << std::endl;
    }
}

void Init_Music(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = ((command >> 24) & 0XFF);
    result = i2c_smbus_write_byte(file, cmd_code);
    if (result == 0)
    {
        std::cout << "Music is initialized successfully." << std::endl;
    }
}

void Load_Music_Config_EPRM(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 cmd_data;
    __s32 result;
    cmd_code = ((command >> 24) & 0XFF);
    cmd_data = ((command >> 16) & 0xFF);
    result = i2c_smbus_write_byte_data(file, cmd_code, cmd_data);
    if (result == 0 && cmd_data == 0x01)
    {
        std::cout << "Bottom Music Configuration is loaded from EEPROM." << std::endl;
    }
    else if (result == 0 && cmd_data == 0x02)
    {
        std::cout << "Top Music Configuration is loaded from EEPROM." << std::endl;
    }
    else
    {
        std::cout << "An error occured. Did you enter the correct music pin?" << std::endl;
    }
}

void Print_Loaded_Music_Config(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 cmd_data;
    __s32 result;
    cmd_code = ((command >> 24) & 0XFF);
    cmd_data = ((command >> 16) & 0xFF);
    result = i2c_smbus_write_byte_data(file, cmd_code, cmd_data);
    if (result == 0 && cmd_data == 0x01)
    {
        std::cout << "Bottom Music Configuration is sent to slave output." << std::endl;
    }
    else if (result == 0 && cmd_data == 0x02)
    {
        std::cout << "Top Music Configuration is sent to slave output." << std::endl;
    }
    else
    {
        std::cout << "An error occured. Did you enter the correct music pin?" << std::endl;
    }
}

void Write_Loaded_Config_to_Music(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 cmd_data;
    __s32 result;
    cmd_code = ((command >> 24) & 0XFF);
    cmd_data = ((command >> 16) & 0xFF);
    result = i2c_smbus_write_byte_data(file, cmd_code, cmd_data);
    if (cmd_data == 0x01)
    {
        std::cout << "Configuration loaded from EEPROM is written to bottom Music." << std::endl;
    }
    else if (cmd_data == 0x02)
    {
        std::cout << "Configuration loaded from EEPROM is written to top Music." << std::endl;
    }
    else
    {
        std::cout << "An error occured. Did you enter the correct music pin?" << std::endl;
    }
}

__u16 Read_from_Music_Register(int file, __u32 command)
{
    __u8 cmd_code;
    __u16 cmd_data;
    __u16 ret_value;
    __s32 result;
    cmd_code = ((command >> 24) & 0XFF);
    cmd_data = ((command >> 8) & 0xFFFF);
    result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    ret_value = result;
    if (((cmd_data >> 8) & 0xFF) == 0x01)
    {
        std::cout << "Bottom Music register is read." << std::endl;
        std::cout << "The result is:" << ret_value << std::endl;
    }
    else if (((cmd_data >> 8) & 0xFF) == 0x02)
    {
        std::cout << "Top Music register is read." << std::endl;
        std::cout << "The result is:" << ret_value << std::endl;
    }
    else
    {
        std::cout << "An error occured. Did you enter the correct music pin?" << std::endl;
    }

    return ret_value;
}

void Write_to_Music_Register(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 cmd_data[3];
    __s32 result;
    cmd_code = ((command >> 24) & 0xFF);
    cmd_data[0] = ((command >> 16) & 0xFF);
    cmd_data[1] = (command & 0xFF);             // This is the LSB of value we are writing to Music
    cmd_data[2] = ((command >> 8) & 0xFF);      // This is the MSB of value we are writing to Music
    result = i2c_smbus_write_block_data(file, cmd_code, 3, cmd_data);
    if ((result == 0) && ((cmd_data[0] >> 6) == 0x01))
    {
        //std::cout << "We wrote to the bottom Music register successfully." << std::endl;
    }
    else if ((result == 0) && ((cmd_data[0] >> 6) == 0x02))
    {
        //std::cout << "We wrote to the top Music register successfully." << std::endl;
    }
    else
    {
        std::cout << "An error occured. Did you enter the correct music pin?" << std::endl;
    }
}

void Set_ADC_Rate(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 cmd_data;
    __s32 result;
    cmd_code = ((command >> 24) & 0XFF);
    cmd_data = ((command >> 16) & 0xFF);
    result = i2c_smbus_write_byte_data(file, cmd_code, cmd_data);
    if (result == 0)
    {
        std::cout << "ADC data rate is updated." << std::endl;
    }
}

void Initialize_ADC(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    if (result == 0)
    {
        std::cout << "ADC is initialized." << std::endl;
    }
}

void Enable_ADC_Run(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    if (result == 0)
    {
        std::cout << "ADC Interrupt Service Routine is activated." << std::endl;
    }
}

void Disable_ADC_Run(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    if (result == 0)
    {
        std::cout << "ADC Interrupt Service Routine is disabled." << std::endl;
    }
}

double Read_Single_Current(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 cmd_data;
    __s32 result;
    double ret_value;
    cmd_code = ((command >> 24) & 0xFF);
    cmd_data = ((command >> 16) & 0xFF);
    result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    ret_value = ((result << 8) | 0x00) * 0.000149011612;
    //std::cout << "ADC Sample from channel " << unsigned(cmd_data) << " is received: " << ret_value << std::endl;
    std::cout << ret_value << std::endl;
    usleep(10000);
    return ret_value;
}

double* Read_All_Current(int file, __u32 command)
{
    __u8 cmd_code;
    __u32 cmd_tmp;
    static double current[16];
    cmd_code = ((command >> 24) & 0xFF);
    for (int i=0; i < 16; i++)
    {
    cmd_tmp = ((CMD_READ_SINGLE_CURRENT << 24) | (i << 16) | 0x0000);
    current[i] = Read_Single_Current(file, cmd_tmp);
    cmd_tmp=0;
    }
    return current;

/*	
//  This part needs to be fixed
    __u8 cmd_code;
    __u8 incoming_bytes[32];
    __u16 current[16];
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_read_block_data(file, cmd_code, incoming_bytes);
    for (int i=0; i < 16; i++)
    {
        current[i] = (((__u16) incoming_bytes[i+1] << 8) | ((__u16) incoming_bytes[i+1]));
    }
    for (int i=0; i < 16; i++)
    {
        std::cout << "Channel " << i << " current: " << current[i] << std::endl;
    }
    std::cout << "ADC Interrupt Service Routine is activated." << std::endl;
    std::cout << result << " bytes have been received." << std::endl;
*/
}

double* Read_ADC_Offset(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 cmd_data;
    __s32 result;
    static double offset[16];
    cmd_code = ((command >> 24) & 0xFF);
    for (__u8 i=0; i < 16; i++)
    {
    cmd_data = i;
    result = i2c_smbus_process_call(file, cmd_code, cmd_data);
    offset[i] = ((result << 8) | 0x00) * 0.000149011612;
    std::cout << "ADC Offset from channel " << unsigned(cmd_data) << " is received: " << offset[i] << std::endl;
    }
    return offset;
}

void Record_ADC_Offset(int file, __u32 command)
{
    bool all_offset_ok = false;
    
    while(!all_offset_ok)
    {   
        all_offset_ok = true;
        __u8 cmd_code;
        __s32 result;
        cmd_code = command >> 24;
        result = i2c_smbus_write_byte(file, cmd_code);
    
        double *rec_offset;
        rec_offset = Read_ADC_Offset(file, 0x1D000000);
        for (int i=0; i<16; i++)
        {
            if(rec_offset[i] < 100.0)
            {
                 all_offset_ok = false;
                 i = 16;
            }
        }
    }
    std::cout << "ADC Offset is recorded correctly." << std::endl;
}

void Turn_on_HV(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "HV Switch is on." << std::endl;    
}

void Turn_off_HV(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "HV Switch is off." << std::endl;    
}

void Enable_HV_CTR(int file, __u32 command)
{
    __u8 cmd_code;
    __s32 result;
    cmd_code = command >> 24;
    result = i2c_smbus_write_byte(file, cmd_code);
    std::cout << "HV control is enabled." << std::endl;    
}

__u8 Check_HV_Status(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 HV_Status;
    __s32 result;
    cmd_code = ((command >> 24) & 0xFF);
    result = i2c_smbus_process_call(file, cmd_code, 0x00);
    HV_Status = (result & 0xFF);
    if(HV_Status == 1) { std::cout << "HV is off." << std::endl;}
    if(HV_Status == 2) { std::cout << "HV is on." << std::endl;}
    return result;
}

void Set_Max_Current(int file, __u32 command)
{
    __u8 cmd_code;
    __u16 cmd_data;
    __s32 result;
    cmd_code = ((command >> 24) & 0xFF);
    cmd_data = ((command >> 8) & 0xFFFF);
    result = i2c_smbus_process_call(file, cmd_code, cmd_data);
}

void Update_Music_Config(int file, __u32 command)
{
    __u8 cmd_code;
    __u8 music_ID;
    __u32 cmd_rebuilt;
    int address[23];
    int value[23];
    int write_cmd = 21;
    std::string filename= "Music_Register_Table.csv";
    Read_Music_Config(filename, address, value);

    cmd_code = ((command >> 24) & 0xFF);
    music_ID = ((command >> 16) & 0xFF);

    for (int i=0; i<23; i++)
    {
        int value_LSB = (value[i] & 0xFF);
        int value_MSB = ((value[i] >> 8) & 0xFF);
        cmd_rebuilt = (((__u32)write_cmd << 24) | ((__u32)music_ID << 22) | ((__u32)address[i] << 16) | ((__u32)value_MSB << 8) | (__u32)(value_LSB));
        Write_to_Music_Register(file, cmd_rebuilt);
        usleep(10000);
    }
    std::cout << "MUSIC Configured" << std::endl;
}

void Read_Music_Config(std::string ConfigCSV, int *address, int *value)
{
    std::ifstream MusicConfig(ConfigCSV);
    if(!MusicConfig.is_open()) throw std::runtime_error("Can not open Music Configuration file.");

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
}

void Set_Bias_Voltage(int i, int file, __u32 command)
{
    std::string filename = "DVOffset.csv";
    std::ifstream BiasVoltage(filename);
    if(!BiasVoltage.is_open()) throw std::runtime_error("Can not open Bias Voltage Configuration file.");

    int startline = i*16;
    int BV_value[16] = {0};
    std::string BV_str;

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
    for (int i=0; i<8; i++)
    {
        cmd_rebuilt = (((__u32)0x15 << 24) | ((__u32)(0x01) << 22) | ((__u32)i << 16) | ((__u32)BV_value[i]));
        std::cout << std::hex << cmd_rebuilt << "\t";
        Write_to_Music_Register(file, cmd_rebuilt);
        usleep(100000);
        cmd_rebuilt = (((__u32)0x15 << 24) | ((__u32)(0x02) << 22) | ((__u32)i << 16) | ((__u32)BV_value[i+8]));
        std::cout << std::hex << cmd_rebuilt << std::endl;
        Write_to_Music_Register(file, cmd_rebuilt);
        usleep(100000);
    }
    std::cout << "Finished Setting bias voltage for SIAB# " << i << std::endl;
}
