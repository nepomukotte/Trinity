#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "anyio.h"
#include "lbp16.h"
#include "eth_boards.h"
#include "TB_Map.h"
#include <mqueue.h>

#define CLIENT_QUEUE_NAME_TB    "/trigger_queue"
#define QUEUE_PERMISSIONS 0660
#define MAX_MESSAGES 5
#define MAX_MSG_SIZE 256
#define MSG_BUFFER_SIZE 256

static board_access_t access_board;

struct Counters_Struct
{
	u64 ClockCounter;
	u32 EventCounter;

	u32 BifocalCounter;
	u32 DiscCounter;
	u32 InternalCounter;
	u32 ExternalCounter;
	u32 GPSCounter;
	u32 HledCounter;

	u32 BifocalLiveCounter;
	u32 DiscLiveCounter;
	u32 InternalLiveCounter;
	u32 ExternalLiveCounter;
	u32 GPSLiveCounter;
	u32 HledLiveCounter;

	u64 TB_BusyCounter;
	u64 CoBo_BusyCounter;
	u64 Transit_BusyCounter;
	u64 Buffer_BusyCounter;
	u64 Mem_FullCounter;
	u64 Mem_ReadingCounter;
	u64 Mem_WritingCounter;

	u32 IO_rate_Counter[64];
	u64 IO_total_counts;
};

int read_32bit(u16 address, u32 data[], u8 size);
int write_32bit(u16 address, u32 data);
void Get_Counters(struct Counters_Struct *Counters);
void Log_Counters(struct Counters_Struct *Counters, int threshold_value);

int main(int argc, char *argv[]) {

    mqd_t qd_server, qd_client;
    struct mq_attr attr;
    char in_buffer [MSG_BUFFER_SIZE];
	struct Counters_Struct Counters;

    attr.mq_flags = 0;
    attr.mq_maxmsg = MAX_MESSAGES;
    attr.mq_msgsize = MAX_MSG_SIZE;
    attr.mq_curmsgs = 0;

    if ((qd_client = mq_open (CLIENT_QUEUE_NAME_TB, O_RDONLY | O_CREAT, QUEUE_PERMISSIONS, &attr)) == -1)
    {
        perror ("Client: mq_open (client)");
        exit (1);
    }

	access_board.eth = 1;
	access_board.address = 1;
	access_board.dev_addr = "192.168.1.121";
	access_board.device_name = "7I80";

	// Initialize, scan and list ETH boards
	eth_boards_init(&access_board);
	eth_boards_scan(&access_board);
	anyio_list_dev(&access_board);

	// Clearing the Busy Bits
	write_32bit(Clear_Busy, All_Bits_High);

	// Enable the Block Memory
	write_32bit(Block_Memory_Enable, All_Bits_High);

	// Set the Block Memory Start Address to zero
	write_32bit(Memory_Start_Address, All_Bits_Low);

	// Reset the Block Memory
	write_32bit(Readout_Done, All_Bits_High);

	// Enabling the Internal Trigger only
	write_32bit(Enable_Trigger_Types, En_Int_Trig_only);

	// Setting the Internal Trigger to 10 Hz
	write_32bit(Internal_Trigger_Prescale, Trig_Rate_10Hz);

	// Setting the Internal Trigger to HLED mode
	write_32bit(Internal_Trigger_Mode, All_Bits_High);

	// Setting the rate counter period
	write_32bit(Rate_Counter_Period, Thirty_Second_Interval);

	// Clearing the Counters
	write_32bit(Clear_Counters, All_Bits_High);

	
	while (1)
	{
		if (mq_receive (qd_client, in_buffer, MSG_BUFFER_SIZE, NULL) == -1)
		{
			perror ("Client: mq_receive");
			exit (1);
		}

		const char *threshold_ptr = in_buffer;
		int threshold_value = (int)strtol(threshold_ptr, NULL, 0);

		if(threshold_value == 0) 
		{
			// Clearing the Counters
			write_32bit(Clear_Counters, All_Bits_High);
			
			printf ("TB: Enabling the trigger\n");
			write_32bit(Enable_Triggering, All_Bits_High);
		}
		else if((threshold_value > 0) && (threshold_value < 512))
		{
			printf ("TB: Disabling the trigger\n");
			write_32bit(Enable_Triggering, Dis_All_Triggers);

			printf ("TB: Saving trigger file with threshold_value = %f\n\n", (threshold_value*1.58));
			write_32bit(Save_Counters, All_Bits_High);
			Get_Counters(&Counters);
			Log_Counters(&Counters, threshold_value);
			write_32bit(Clear_Counters, All_Bits_High);
			write_32bit(Readout_Done, All_Bits_High);
		}
	}

	return 0;
}

int read_32bit(u16 address, u32 data[], u8 size)
{
	int send = 0; int recv = 0;
	lbp16_cmd_addr packet;
	u32 command = ((0x4200) | (0x80 | size));
	LBP16_INIT_PACKET4(packet, command, address);
	send = eth_socket_send_packet(&packet, sizeof(packet));
	sleep_ns(5*1000*1000);
	recv = eth_socket_recv_packet(data, size*sizeof(data[0]));
	return recv;
}

int write_32bit(u16 address, u32 data)
{
	int send = 0;
	lbp16_cmd_addr_data32 packet;
	LBP16_INIT_PACKET8(packet, 0xC281, address, data);
	send = eth_socket_send_packet(&packet, sizeof(packet));
	sleep_ns(5*1000*1000);
	return send;
}

void Get_Counters(struct Counters_Struct *Counters)
{
	// Reading 29 registers starting from Counters_First_Addr (0x2020)
	// This array will be taken apart to calculate individual counters
	u32 tmp[29];
	read_32bit(Counters_First_Addr, tmp, 29);

	Counters->ClockCounter = (((u64)tmp[1] << 32 )| (u64)tmp[0]);
	Counters->EventCounter = tmp[2];

	Counters->BifocalCounter 	= tmp[3];
	Counters->DiscCounter 	 	= tmp[4];
	Counters->InternalCounter 	= tmp[5];
	Counters->ExternalCounter 	= tmp[6];
	Counters->GPSCounter 	 	= tmp[7];
	Counters->HledCounter 	 	= tmp[8];

	Counters->BifocalLiveCounter	= tmp[9];
	Counters->DiscLiveCounter		= tmp[10];
	Counters->InternalLiveCounter	= tmp[11];
	Counters->ExternalLiveCounter	= tmp[12];
	Counters->GPSLiveCounter		= tmp[13];
	Counters->HledLiveCounter		= tmp[14];

	Counters->TB_BusyCounter		= (((u64)tmp[16] << 32 )| (u64)tmp[15]);
	Counters->CoBo_BusyCounter		= (((u64)tmp[18] << 32 )| (u64)tmp[17]);
	Counters->Transit_BusyCounter	= (((u64)tmp[20] << 32 )| (u64)tmp[19]);
	Counters->Buffer_BusyCounter	= (((u64)tmp[22] << 32 )| (u64)tmp[21]);
	Counters->Mem_FullCounter		= (((u64)tmp[24] << 32 )| (u64)tmp[23]);
	Counters->Mem_ReadingCounter	= (((u64)tmp[26] << 32 )| (u64)tmp[25]);
	Counters->Mem_WritingCounter	= (((u64)tmp[28] << 32 )| (u64)tmp[27]);

	// Reading 64 registers starting from IO_0_rate_Addr (0x2094)
	// This array will be taken apart to calculate individual IO rate counters
	u32 tmp_IO[64];
	read_32bit(IO_0_rate_Addr, tmp_IO, 64);
	Counters->IO_total_counts = 0;
	for (int i=0; i<64; i++)
	{
		Counters->IO_rate_Counter[i] = tmp_IO[i];
		Counters->IO_total_counts += tmp_IO[i];
	}
}

void Log_Counters(struct Counters_Struct *Counters, int threshold_value)
{
	printf("Here are the counters value:\n");

	printf("Clock: %lu\n", Counters->ClockCounter);
	printf("Event: %u\n", Counters->EventCounter);

	printf("Bifocal: %d\n", Counters->BifocalCounter);
	printf("Disc Test: %d\n", Counters->DiscCounter);
	printf("Internal: %d\n", Counters->InternalCounter);
	printf("External: %d\n", Counters->ExternalCounter);
	printf("GPS: %d\n", Counters->GPSCounter);
	printf("HLED: %d\n", Counters->HledCounter);
	
	printf("Bifocal Live: %d\n", Counters->BifocalLiveCounter);
	printf("Disc Test Live: %d\n", Counters->DiscLiveCounter);
	printf("Internal Live: %d\n", Counters->InternalLiveCounter);
	printf("External Live: %d\n", Counters->ExternalLiveCounter);
	printf("GPS Live: %d\n", Counters->GPSLiveCounter);
	printf("HLED Live: %d\n", Counters->HledLiveCounter);

	printf("TB Busy: %lu\n", Counters->TB_BusyCounter);
	printf("CoBo Busy: %lu\n", Counters->CoBo_BusyCounter);
	printf("Transit Busy: %lu\n", Counters->Transit_BusyCounter);
	printf("Buffer Busy: %lu\n", Counters->Buffer_BusyCounter);
	printf("Memory Full: %lu\n", Counters->Mem_FullCounter);
	printf("Memory Reading: %lu\n", Counters->Mem_ReadingCounter);
	printf("Memory Writing: %lu\n", Counters->Mem_WritingCounter);

	FILE *fileText;
	fileText = fopen("./IO_Counts.txt","a+");

	fprintf(fileText, "%f\n",threshold_value*1.58);
	for (int i=0; i<64; i++)
	{
		printf("IO#%d: %d\n",i, Counters->IO_rate_Counter[i]);
		fprintf(fileText, "%d %d\n",i, Counters->IO_rate_Counter[i]);
	}
	fprintf(fileText, "\n");
	fclose(fileText);
	printf("Total IO Counts: %ld\n",Counters->IO_total_counts);
}
