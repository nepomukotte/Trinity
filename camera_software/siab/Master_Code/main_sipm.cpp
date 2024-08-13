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
	int total_files = 24;
	int total_samples = 60;
	time_t seconds;

	__u32 Packet_READ_SiPM_Temp = 0x0A000000;

	file = open("/dev/i2c-1", O_RDWR);
	if (file < 0) {
		std::cout << "Error in opening file ..." << std::endl;
		exit(1);
	}

	while(1) {
		for(int i=0; i<total_files; i++) {
			FILE *SiPMfilePointer;
			char SiPMfilename[20];
			sprintf(SiPMfilename, "file_SiPM_%d.txt", i);
			SiPMfilePointer = fopen(SiPMfilename, "wb");
			for(int j=0; j<total_samples; j++){
				__u16 sipm_temp[NofSIABs];

				for (int k=0; k<NofSIABs; k++) {
					I2C_ADDR = First_SIAB_ID + k;
					if (ioctl(file, I2C_SLAVE, I2C_ADDR) < 0) {
						std::cout << "Error in setting up SIAB Number: " << k+1 << std::endl;
						exit(1);
					}
					sipm_temp[k] = Read_SiPM_Temp(file, Packet_READ_SiPM_Temp);
					usleep(10000);
				}

				seconds = time (NULL);
				std::cout << "seconds since January 1, 1970: " << seconds << std::endl;
				fwrite(&seconds, sizeof(__u16), 1, filePointer);

				fwrite(sipm_temp, sizeof(__u16), NofSIABs, SiPMfilePointer);
				sleep(5);
			}
			fclose(SiPMfilePointer);
			std::cout << "Finished SiPM File# " << i+1 << std::endl;
		}
	}

	return 0;
}
