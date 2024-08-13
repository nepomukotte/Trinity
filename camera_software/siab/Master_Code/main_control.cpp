#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/i2c-dev.h>
#include <linux/i2c.h>
#include <linux/smbus.h>
#include "commands.h"

int main (void)
{
	int file = 0;
	int I2C_ADDR = 0;
	int First_SIAB_ID = 0x10;
	double* Current_all;
	unsigned int temperature;
	std::string user_input;

	file = open("/dev/i2c-1", O_RDWR);
	if (file < 0)
	{
		std::cout << "Error in opening file ..." << std::endl;
		exit(1);
	}

	while ( std::getline(std::cin, user_input))
	{
		__u8 CMD_ID;
		__u32 SIAB_ID;
		__u32 CMD_Packet;
		__u64 Packet;

		std::stringstream user_string;
		user_string << user_input;
		user_string >> std::hex >> Packet;

		SIAB_ID = ((Packet >> 32) & 0xFFFFFFFF);
		CMD_ID = ((Packet >> 24) & 0xFF);
		CMD_Packet = (Packet & 0xFFFFFFFF);

		for (int i=0; i<32; i++)
		{
			int SIAB_bit = 0;
			SIAB_bit = ((SIAB_ID >> i) & 0x01);
			if (SIAB_bit)
			{
				I2C_ADDR = First_SIAB_ID + i;
				if (ioctl(file, I2C_SLAVE, I2C_ADDR) < 0)
				{
					std::cout << "Error in setting up SIAB Number: " << i+1 << std::endl;
					exit(1);
				}

				switch (CMD_ID)
				{
					case CMD_RESET:
						Reset_Atmega328P(file, CMD_Packet);
						break;
					case CMD_TURN_ON_LED:
						Turn_on_LED(file, CMD_Packet);
						break;
					case CMD_TURN_OFF_LED:
						Turn_off_LED(file, CMD_Packet);
						break;
					case CMD_TURN_ON_3V3:
						Turn_on_3V3(file, CMD_Packet);
						break;
					case CMD_TURN_OFF_3V3:
						Turn_off_3V3(file, CMD_Packet);
						break;
					case CMD_TURN_ON_5V:
						Turn_on_5V(file, CMD_Packet);
						break;
					case CMD_TURN_OFF_5V:
						Turn_off_5V(file, CMD_Packet);
						break;
					case CMD_TURN_ON_3V3_5V:
						Turn_on_3V3_5V(file, CMD_Packet);
						break;
					case CMD_TURN_OFF_3V3_5V:
						Turn_off_3V3_5V(file, CMD_Packet);
						break;
					case CMD_READ_SIPM_TEMP:
						Read_SiPM_Temp(file, CMD_Packet);
						break;
					case CMD_READ_UC_TEMP:
						Read_UC_Temp(file, CMD_Packet);
						break;
					case CMD_READ_EPRM_BYTE:
						Read_EPRM_Byte(file, CMD_Packet);
						break;
					case CMD_READ_EPRM_WORD:
						Read_EPRM_Word(file, CMD_Packet);
						break;
					case CMD_WRITE_EPRM_BYTE:
						Write_EPRM_Byte(file, CMD_Packet);
						break;
					case CMD_WRITE_EPRM_WORD:
						Write_EPRM_Word(file, CMD_Packet);
						break;
					case CMD_INIT_MUSIC:
						Init_Music(file, CMD_Packet);
						break;
					case CMD_LOAD_MUSIC_CONFIG_EPRM:
						Load_Music_Config_EPRM(file, CMD_Packet);
						break;
					case CMD_PRINT_LOADED_MUSIC_CONFIG:
						Print_Loaded_Music_Config(file, CMD_Packet);
						break;
					case CMD_WRITE_LOADED_CONFIG_TO_MUSIC:
						Write_Loaded_Config_to_Music(file, CMD_Packet);
						break;
					case CMD_READ_FROM_MUSIC_REGISTER:
						Read_from_Music_Register(file, CMD_Packet);
						break;
					case CMD_WRITE_TO_MUSIC_REGISTER:
						Write_to_Music_Register(file, CMD_Packet);
						break;
					case CMD_SET_ADC_RATE:
						Set_ADC_Rate(file, CMD_Packet);
						break;
					case CMD_INITIALIZE_ADC:
						Initialize_ADC(file, CMD_Packet);
						break;
					case CMD_ENABLE_ADC_RUN:
						Enable_ADC_Run(file, CMD_Packet);
						break;
					case CMD_DISABLE_ADC_RUN:
						Disable_ADC_Run(file, CMD_Packet);
						break;
					case CMD_READ_SINGLE_CURRENT:
						Read_Single_Current(file, CMD_Packet);
						break;
					case CMD_READ_ALL_CURRENT:
						Current_all = Read_All_Current(file, CMD_Packet);
						break;
					case CMD_RECORD_ADC_OFFSET:
						Record_ADC_Offset(file, CMD_Packet);
						break;
					case CMD_READ_ADC_OFFSET:
						Read_ADC_Offset(file, CMD_Packet);
						break;
					case CMD_TURN_ON_HV:
						Turn_on_HV(file, CMD_Packet);
						break;
					case CMD_TURN_OFF_HV:
						Turn_off_HV(file, CMD_Packet);
						break;
					case CMD_ENABLE_HV_CTR:
						Enable_HV_CTR(file, CMD_Packet);
						break;
					case CMD_CHECK_HV_STATUS:
						Check_HV_Status(file, CMD_Packet);
						break;
					case CMD_SET_MAX_CURRENT:
						Set_Max_Current(file, CMD_Packet);
						break;
					case CMD_UPDATE_MUSIC_CONFIG:
						Update_Music_Config(file, CMD_Packet);
						break;
					case CMD_SET_BIAS_VOLTAGE:
						Set_Bias_Voltage(i, file, CMD_Packet);
						break;
					default:
						break;
				}
			}
		}
	}

	return 0;
}
