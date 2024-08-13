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

#define CLIENT_QUEUE_NAME_TB	"/trigger_queue"
#define MSG_BUFFER_SIZE 256

int main()
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

	mqd_t qd_client;
	std::string user_input;
	char in_buffer [MSG_BUFFER_SIZE];
	std::string start_buffer = "0";
	char zero_buffer[2];
	strcpy(zero_buffer, start_buffer.c_str());

	file = open("/dev/i2c-0", O_RDWR);
	if (file < 0)
	{
		std::cout << "Error in opening file ..." << std::endl;
		exit(1);
	}

	if ((qd_client = mq_open (CLIENT_QUEUE_NAME_TB, O_WRONLY)) == -1)
	{
		perror ("Server: Not able to open client queue");
		exit (1);
	}

	for (int i=0; i<NofSteps; i++)
	{
		threshold_value = threshold_start + i*threshold_step_size;
		int threshold_reg_value = ((threshold_value << 3) | 0xF007);
		std::cout << "Controller: Starting to set the discriminator threshold at: " << (threshold_value*1.58) << " mV" << std::endl;

		std::string buffer;
		buffer = std::to_string(threshold_value);
		int n = buffer.length();
		char out_buffer[n+1];
		strcpy(out_buffer, buffer.c_str());

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
			for (int i=8; i<NofChannels+8; i++)
			{
				CMD_Packet = ((0x15 << 24) | (((0x01 << 6) | i) << 16) | threshold_reg_value);
				Write_to_Music_Register(file, CMD_Packet);
				usleep(10000);
				CMD_Packet = ((0x15 << 24) | (((0x02 << 6) | i) << 16) | threshold_reg_value);
				Write_to_Music_Register(file, CMD_Packet);
			}
		}
		std::cout << "Controller: Finished setting the discriminator threshold for all SIABs." << std::endl;

		if (mq_send (qd_client, zero_buffer, strlen (zero_buffer) + 1, 0) == -1)
		{
			perror ("Unable to send start message to trigger board");
			continue;
		}
		cout << "Controller: Enable trigger message was sent to trigger board." << std::endl;

		sleep(5);

		if (mq_send (qd_client, out_buffer, strlen (out_buffer) + 1, 0) == -1)
		{
			perror ("Unable to send stop message to trigger board");
			continue;
		}
		cout << "Controller: Disable trigger message was sent to trigger board.\n" << std::endl;
	}

	return 0;
}