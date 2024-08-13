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

int main (void)
{
	int file = 0;
	int I2C_ADDR = 0;
	int NofSIABs = 16;
	int First_SIAB_ID = 0x10;
	time_t seconds;

	__u32 Packet_READ_Current = 0x1A090000;
	__u32 Packet_Check_HV = 0x21000000;

	file = open("/dev/i2c-1", O_RDWR);
	if (file < 0) {
		std::cout << "Error in opening file ..." << std::endl;
		exit(1);
	}

	while(1)
	{
		I2C_ADDR = First_SIAB_ID + 8;
		if (ioctl(file, I2C_SLAVE, I2C_ADDR) < 0)
		{
			std::cout << "Error in setting up SIAB Number 1" << std::endl;
			exit(1);
		}
		Read_Single_Current(file, Packet_READ_Current);
		Check_HV_Status(file, Packet_Check_HV);
		std::cout << std::endl;
		sleep(2);
	}
	return 0;
}
