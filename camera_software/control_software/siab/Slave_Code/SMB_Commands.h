#ifndef __SMB_COMMANDS_H__
#define __SMB_COMMANDS_H__

#include "SMB_Slave.h"

const int Music_3V3_Enable = 17;    // PC3
const int Music_5V_Enable = 16;     // PC2
const int LED_Pin = 2;              // PD2
const int HV_Enable = 3;            // PD3
const int Temp_Pin = A7;            // ADC7

#define CMD_RESET                 0x01
#define CMD_TURN_ON_LED           0x02
#define CMD_TURN_OFF_LED          0x03
#define CMD_TURN_ON_3V3           0x04
#define CMD_TURN_OFF_3V3          0x05
#define CMD_TURN_ON_5V            0x06
#define CMD_TURN_OFF_5V           0x07
#define CMD_TURN_ON_3V3_5V        0x08
#define CMD_TURN_OFF_3V3_5V       0x09

#define CMD_READ_SIPM_TEMP        0x0A
#define CMD_READ_UC_TEMP          0X0B

#define CMD_READ_EPRM_BYTE        0X0C
#define CMD_READ_EPRM_WORD        0X0D
#define CMD_WRITE_EPRM_BYTE       0X0E
#define CMD_WRITE_EPRM_WORD       0X0F

#define CMD_INIT_MUSIC                      0x10
#define CMD_LOAD_MUSIC_CONFIG_EPRM          0x11
#define CMD_PRINT_LOADED_MUSIC_CONFIG       0x12
#define CMD_WRITE_LOADED_CONFIG_TO_MUSIC    0x13
#define CMD_READ_FROM_MUSIC_REGISTER        0x14
#define CMD_WRITE_TO_MUSIC_REGISTER         0x15

#define CMD_SET_ADC_RATE          0x16
#define CMD_INITIALIZE_ADC        0x17
#define CMD_ENABLE_ADC_RUN        0x18
#define CMD_DISABLE_ADC_RUN       0x19
#define CMD_READ_SINGLE_CURRENT   0x1A
#define CMD_READ_ALL_CURRENT      0x1B
#define CMD_RECORD_ADC_OFFSET     0x1C
#define CMD_READ_ADC_OFFSET       0x1D

#define CMD_TURN_ON_HV            0x1E
#define CMD_TURN_OFF_HV           0x1F
#define CMD_ENABLE_HV_CTR         0x20
#define CMD_CHECK_HV_STATUS       0x21
#define CMD_SET_MAX_CURRENT       0x22

#define CMD_GET_SIAB_ID           0x25

void ProcessReceiveByte(SMBData *smb);
void ProcessMessage(SMBData *smb);

void Reset_Atmega328P(SMBData *smb);
void Turn_on_LED(SMBData *smb);
void Turn_off_LED(SMBData *smb);
void Turn_on_3V3(SMBData *smb);
void Turn_off_3V3(SMBData *smb);
void Turn_on_5V(SMBData *smb);
void Turn_off_5V(SMBData *smb);
void Turn_on_3V3_5V(SMBData *smb);
void Turn_off_3V3_5V(SMBData *smb);

void Read_SiPM_Temp(SMBData *smb);
void Read_UC_Temp(SMBData *smb);

void Read_EPRM_Byte(SMBData *smb);
void Read_EPRM_Word(SMBData *smb);
void Write_EPRM_Byte(SMBData *smb);
void Write_EPRM_Word(SMBData *smb);

void Init_Music(SMBData *smb);
void Load_Music_Config_EPRM(SMBData *smb);
void Print_Loaded_Music_Config(SMBData *smb);
void Write_Loaded_Config_to_Music(SMBData *smb);
void Read_from_Music_Register(SMBData *smb);
void Write_to_Music_Register(SMBData *smb);

void Set_ADC_Rate(SMBData *smb);
void Initialize_ADC(SMBData *smb);
void Enable_ADC_Run(SMBData *smb);
void Disable_ADC_Run(SMBData *smb);
void Read_Single_Current(SMBData *smb);
void Read_All_Current(SMBData *smb);
void Record_ADC_Offset(SMBData *smb);
void Read_ADC_Offset(SMBData *smb);

void Turn_on_HV(SMBData *smb);
void Turn_off_HV(SMBData *smb);
void Enable_HV_CTR(SMBData *smb);
void Check_HV_Status(SMBData *smb);
void Set_Max_Current(SMBData *smb);

void Get_SIAB_ID(SMBData *smb);
void ADC_ready(const char* message, bool pinstate);

#endif
