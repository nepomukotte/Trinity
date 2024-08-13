#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <time.h>
#include <sys/ioctl.h>
#include <linux/i2c-dev.h>
#include <linux/i2c.h>
#include <linux/smbus.h>
#include "commands.h"

using namespace std;

int main (int argc, char *argv[])
{
	int file = 0;
	int I2C_ADDR = 0;
	int NofSIABs = 32;
	int SIAB_ID_List[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,32,31};
	int length = sizeof(SIAB_ID_List)/sizeof(SIAB_ID_List[0]);
	int First_SIAB_ID = 0x10;
	uint32_t select_siab = 0;
	std::string selected_siab_str = argv[1];
	std::stringstream user_input;
	user_input << selected_siab_str;
	user_input >> std::hex >> select_siab;

	__u32 Packet_Turn_On_Music  	= 0x08000000;
	__u32 Packet_Init_Music     	= 0x10000000;
	__u32 Packet_Config_Bot_Music   = 0x23010000;
	__u32 Packet_Config_Top_Music   = 0x23020000;

	file = open("/dev/i2c-0", O_RDWR);
	if (file < 0) {
		std::cout << "Error in opening file ..." << std::endl;
		exit(1);
	}
	
	for (int k=0; k<NofSIABs; k++)
	{
		int SIAB_bit = 0;
		SIAB_bit = ((select_siab >> k) & 0x01);
		if(SIAB_bit)
		{
			std::cout << "SIAB# " << k << std::endl;
			I2C_ADDR = First_SIAB_ID + SIAB_ID_List[k];
			if (ioctl(file, I2C_SLAVE, I2C_ADDR) < 0)
			{
				std::cout << "Error in setting up SIAB Number: " << k+1 << std::endl;
				exit(1);
			}

			usleep(100000);
			Turn_on_3V3_5V(file, Packet_Turn_On_Music);
			usleep(1000000);
			Init_Music(file, Packet_Init_Music);
			usleep(500000);
			Update_Music_Config(file, Packet_Config_Bot_Music);
			usleep(500000);
			Update_Music_Config(file, Packet_Config_Top_Music);		
			usleep(500000);
		}
	}

	return 0;
}