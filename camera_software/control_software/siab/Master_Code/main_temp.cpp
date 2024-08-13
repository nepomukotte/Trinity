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
	int First_SIAB_ID = 0x10;
	int total_files = 1;
	int total_samples = 300000;
	int chamber_temp = atoi(argv[1]);
	int SIAB_ID_List[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,32,31};
	int length = sizeof(SIAB_ID_List)/sizeof(SIAB_ID_List[0]);
	time_t start_time;
	time_t sample_time;

	__u32 Packet_READ_SiPM_Temp = 0x0A000000;
	__u32 Packet_READ_UC_Temp   = 0x0B000000;

	file = open("/dev/i2c-0", O_RDWR);
	if (file < 0) {
		std::cout << "Error in opening file ..." << std::endl;
		exit(1);
	}

	start_time = time (NULL);
	std::cout << std::endl << "Run Started at: " << start_time << std::endl;
	
	for(int i=0; i<total_files; i++) {
		FILE *filePointer;
		char filename[30];
		sprintf(filename, "TVAC/T_Calib_%d_C.txt", chamber_temp);
		filePointer = fopen(filename, "wb");
		for(int j=0; j<total_samples; j++){
			__u32 uc_temp[NofSIABs];
			__u32 sipm_temp[NofSIABs];
			sample_time = time (NULL);
			std::cout << std::endl << "Seconds since start of the run: " << sample_time - start_time << std::endl;
			for (int k=0; k<length; k++) {
				std::cout << "SIAB ID: " << k << "\t";
				I2C_ADDR = First_SIAB_ID + SIAB_ID_List[k];
				if (ioctl(file, I2C_SLAVE, I2C_ADDR) < 0) {
					std::cout << "Error in setting up SIAB ID: " << k << std::endl;
					exit(1);
				}
				usleep(50000);
				uc_temp[k] = Read_UC_Temp(file, Packet_READ_UC_Temp);
				usleep(50000);
				sipm_temp[k] = Read_SiPM_Temp(file, Packet_READ_SiPM_Temp);
				usleep(50000);
			}

			fwrite(&sample_time, sizeof(__u32), 1, filePointer);
			fwrite(uc_temp, sizeof(__u32), NofSIABs, filePointer);
			fwrite(sipm_temp, sizeof(__u32), NofSIABs, filePointer);
			sleep(3);
		}
		fclose(filePointer);
		std::cout << "Finished." << std::endl;
	}

	return 0;
}