#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fcntl.h>
#include <sstream>
#include <unistd.h>
#include <thread>
#include <sys/types.h>
#include <sys/stat.h>
#include <mqueue.h>

#include <linux/i2c-dev.h>
#include <linux/i2c.h>
#include <linux/smbus.h>
#include "commands.h"

using namespace std;

int main(int argc, char** argv)
{
	int file = 0;
	int I2C_ADDR = 0;
	int NofSIAB = 32;
	int NofChannels = 8;
	int First_SIAB_ID = 0x10;
	time_t seconds;

	int NofSteps = 100;
	int threshold_value = 0;
	int threshold_start = 1;
	int threshold_step_size = 5;

	cout<<"Trigger DVth: "<<atoi(argv[1])<<endl;

	if(atoi(argv[2])==32){
		cout<<"Siabs: "<<"ALL"<<endl;
	}
	else{
		cout<<"Siab: "<<atoi(argv[2])<<endl;
	}
	
	if(atoi(argv[3])==1){
		cout<<"BOTTOM MUSIC"<<endl;

	}
	else if (atoi(argv[3])==2){
		cout<<"TOP MUSIC"<<endl;
	}
	else{
		cout<<"BOTH MUSICs"<<endl;
	}

	file = open("/dev/i2c-0", O_RDWR);
	if (file < 0)
	{
		std::cout << "Error in opening file ..." << std::endl;
		exit(1);
	}
	
	threshold_value = atoi(argv[1]);
	//cout<<"char"<<endl;
	int threshold_reg_value = ((threshold_value << 3) | 0xF007);
	std::cout << "Controller: Starting to set the discriminator threshold at: " << (threshold_value) << std::endl;
	std::string buffer;
	buffer = std::to_string(threshold_value);
	int n = buffer.length();
	char out_buffer[n+1];
	strcpy(out_buffer, buffer.c_str());

	//cout<<"Hi"<<endl;
	if(atoi(argv[2])>31){
		for (int i=0; i<NofSIAB; i++)
		{
			usleep(100000);
			__u32 CMD_Packet = 0;
			I2C_ADDR = First_SIAB_ID + i;
			if (ioctl(file, I2C_SLAVE, I2C_ADDR) < 0)
			{
				std::cout << "Error in setting up SIAB Number: " << i+1 << std::endl;
				exit(1);
			}
			for (int j=8; j<NofChannels+8; j++)
			{
				CMD_Packet = ((0x15 << 24) | (((0x02 << 6) | j) << 16) | threshold_reg_value);
				//cout<<CMD_Packet<<endl;
				Write_to_Music_Register(file, CMD_Packet);
				usleep(10000);
				CMD_Packet = 0;
				CMD_Packet = ((0x15 << 24) | (((0x01 << 6) | j) << 16) | threshold_reg_value);
				Write_to_Music_Register(file, CMD_Packet);
				usleep(10000);
			}
		}
	}
	else{
		__u32 CMD_Packet = 0;
		I2C_ADDR = First_SIAB_ID + atoi(argv[2]);
		if (ioctl(file, I2C_SLAVE, I2C_ADDR) < 0)
		{
			std::cout << "Error in setting up SIAB Number: " << argv[2] << std::endl;
			exit(1);
		}
		for (int j=8; j<NofChannels+8; j++)
		{	
			if(atoi(argv[3])==1){
				CMD_Packet = ((0x15 << 24) | (((0x01 << 6) | j) << 16) | threshold_reg_value);
				Write_to_Music_Register(file, CMD_Packet);
			}
			else if(atoi(argv[3])==2){
				CMD_Packet = ((0x15 << 24) | (((0x02 << 6) | j) << 16) | threshold_reg_value);
				Write_to_Music_Register(file, CMD_Packet);
			}
			else{
				CMD_Packet = ((0x15 << 24) | (((0x01 << 6) | j) << 16) | threshold_reg_value);
				Write_to_Music_Register(file, CMD_Packet);
				usleep(10000);
				CMD_Packet = ((0x15 << 24) | (((0x02 << 6) | j) << 16) | threshold_reg_value);
				Write_to_Music_Register(file, CMD_Packet);

			}
			
			usleep(10000);
			
		}

	}
	
	std::cout << "Controller: Finished setting the discriminator threshold for all SIABs." << std::endl;

	return atoi(argv[1]);
}
