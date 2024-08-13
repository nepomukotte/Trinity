#ifndef _COMMANDS_H
#define _COMMANDS_H

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <fstream>
#include <fcntl.h>
#include <string.h>
#include <sys/ioctl.h>
#include <linux/i2c-dev.h>
#include <linux/i2c.h>
#include <i2c/smbus.h>

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
#define CMD_SET_MAX_CURRENT		  0x22

#define CMD_UPDATE_MUSIC_CONFIG	  0x23
#define CMD_SET_BIAS_VOLTAGE	  0x24

void Reset_Atmega328P(int file, __u32 command);
void Turn_on_LED(int file, __u32 command);
void Turn_off_LED(int file, __u32 command);
void Turn_on_3V3(int file, __u32 command);
void Turn_off_3V3(int file, __u32 command);
void Turn_on_5V(int file, __u32 command);
void Turn_off_5V(int file, __u32 command);
void Turn_on_3V3_5V(int file, __u32 command);
void Turn_off_3V3_5V(int file, __u32 command);
__u16 Read_SiPM_Temp(int file, __u32 command);
__u16 Read_UC_Temp(int file, __u32 command);

__u8  Read_EPRM_Byte(int file, __u32 command);
__u16 Read_EPRM_Word(int file, __u32 command);
void  Write_EPRM_Byte(int file, __u32 command);
void  Write_EPRM_Word(int file, __u32 command);

void  Init_Music(int file, __u32 command);
void  Load_Music_Config_EPRM(int file, __u32 command);
void  Print_Loaded_Music_Config(int file, __u32 command);
void  Write_Loaded_Config_to_Music(int file, __u32 command);
__u16 Read_from_Music_Register(int file, __u32 command);
void  Write_to_Music_Register(int file, __u32 command);

void Set_ADC_Rate(int file, __u32 command);
void Initialize_ADC(int file, __u32 command);
void Enable_ADC_Run(int file, __u32 command);
void Disable_ADC_Run(int file, __u32 command);
double Read_Single_Current(int file, __u32 command);
double* Read_All_Current(int file, __u32 command);
void Record_ADC_Offset(int file, __u32 command);
double* Read_ADC_Offset(int file, __u32 command);

void Turn_on_HV(int file, __u32 command);
void Turn_off_HV(int file, __u32 command);
void Enable_HV_CTR(int file, __u32 command);
__u8 Check_HV_Status(int file, __u32 command);
void Set_Max_Current(int file, __u32 command);

void Update_Music_Config(int file, __u32 command);
void Read_Music_Config(std::string ConfigCSV, int *address, int *value);
void Set_Bias_Voltage(int i, int file, __u32 command);

#endif
