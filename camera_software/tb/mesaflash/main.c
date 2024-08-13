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

static board_access_t access_board;
static char bitfile_name[255];
static int fallback_flag = 0;

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

void Get_Counters(struct Counters_Struct *Counters);
void Log_Counters(struct Counters_Struct *Counters);

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


int main(int argc, char *argv[]) {

	int ret = 0;
	int program_flag = 0;
	char bitfile_str[] = "TB_VHDL_Code_V1.0.bit";
	strncpy(bitfile_name, bitfile_str, sizeof(bitfile_name));

	u32 NofEvents[1] = {0};
	u32 NofWords = 0;
	int BlockSize = 100;
	struct Counters_Struct Counters;

	FILE *filePointer;
	board_t *board = NULL;

	access_board.eth = 1;
	access_board.address = 1;
	access_board.dev_addr = "192.168.1.121";
	access_board.device_name = "7I80";

	// Initializing ETH boards
	eth_boards_init(&access_board);
	// Scanning for available devices
	eth_boards_scan(&access_board);
	// List all ethernet devices
	anyio_list_dev(&access_board);

	if(program_flag)
	{
		board = anyio_get_dev(&access_board);
		if(board == NULL)
		{
			printf("No board was found\n");
			return 0;
		}
		board->open(board);
		ret = anyio_dev_write_flash(board, bitfile_name, fallback_flag);
		board->close(board);
		anyio_cleanup(&access_board);
		return ret;
	}
	else
	{
	// Clearing the Busy Bits
	write_32bit(Clear_Busy, All_Bits_High);

	// Enable the Block Memory
	write_32bit(Block_Memory_Enable, All_Bits_High);

	// Set the Block Memory Start Address to zero
	write_32bit(Memory_Start_Address, All_Bits_Low);

	// Reset the Block Memory
	write_32bit(Readout_Done, All_Bits_High);

	// Setting the Transit Busy time to 10 clock ticks (100ns)
	write_32bit(Transit_Busy_Length, Hundred_ns_Interval);

// ********** Reading the Counters *************************** //

	// Setting the rate counter period
	write_32bit(Rate_Counter_Period, One_Second_Interval);

	// Clearing the Counters
	write_32bit(Clear_Counters, All_Bits_High);

	sleep(3);

	Get_Counters(&Counters);
	Log_Counters(&Counters);

// ********************************************************** //

/*
// ********** Setting up Test Trigger *********************** //

	// Enabling the Strech logic for all 64 IO
	write_32bit(Disc_Stretch_Enable_0, All_Bits_High);
	write_32bit(Disc_Stretch_Enable_1, All_Bits_High);

	// Setting the Strech length to MAX value
	write_32bit(Disc_Stretch_Length, All_Bits_High);

	// Setting the Test Trigger Pre-scale value to 1
	write_32bit(Disc_Test_Trigger_Prescale, 1);

	// Enabling the Test Trigger only
	write_32bit(Enable_Trigger_Types, En_Disc_Test_Trig_only);

// ********************************************************** //
*/
/*	
// ********** Setting up Internal Trigger ******************* //

	// Enabling the Internal Trigger only
	//write_32bit(Enable_Trigger_Types, En_Int_Trig_only);

	// Setting the Internal Trigger to 10 Hz
	//write_32bit(Internal_Trigger_Prescale, Trig_Rate_100Hz);

	// Setting the Internal Trigger to Standalone mode
	//write_32bit(Internal_Trigger_Mode, All_Bits_Low);

// ********************************************************** //
*/
/*	
// ********** Setting up External Trigger ******************* //

	// Enabling the External Trigger only
	//write_32bit(Enable_Trigger_Types, En_Ext_Trig_only);

	// Setting the External Trigger to Standalone mode
	//write_32bit(External_Trigger_Mode, All_Bits_Low);

// ********************************************************** //
*/
/*	
// ********** Setting up Bi-Focal Trigger ******************* //

	// Enabling the Strech logic for all 64 IO
	write_32bit(Disc_Stretch_Enable_0, All_Bits_High);
	write_32bit(Disc_Stretch_Enable_1, All_Bits_High);

	// Setting the strech length to MAX value
	write_32bit(Disc_Stretch_Length, All_Bits_High);

	// Enabling the Bi-Focal Trigger only
	write_32bit(Enable_Trigger_Types, En_BiFocal_Trig_only);

// ********************************************************** //
*/
	write_32bit(Enable_Triggering, All_Bits_High);

	sleep(3);

	// Disabling the overall Trigger
	write_32bit(Enable_Triggering, All_Bits_Low);

	// Starting the Readout Process
	write_32bit(Readout_Start, 0x00000001);
	sleep_ns(5*1000*1000);

	// Reading the number of events in the Memory
	read_32bit(Events_Written, NofEvents, 1);
	printf("Number of events: %d\n", NofEvents[0]);
	NofWords = (5 * NofEvents[0]);
	printf("Number of Words: %d\n", NofWords);

	// Writing Events_Writtents to a file
	filePointer = fopen("data.bin", "wb");

	if ((NofWords > 0) && (NofWords <= 8192))
	{
		//First piece of Memory with Block Select Regiter = 0
		int NofBlocks = (NofWords / BlockSize);
		int NofLastBlock = (NofWords % BlockSize);
		write_32bit(Memory_Block_Select, All_Bits_Low);
		for (int i=0; i<NofBlocks; i++)
		{
			u32 data_word[BlockSize];
			read_32bit(Block_Memory_Access_Reg_Start_Addr+4*BlockSize*i, data_word, BlockSize);
			fwrite(data_word, sizeof(u32), BlockSize, filePointer);
		}
		if(NofLastBlock != 0)
		{
			u32 last_data_block[NofLastBlock];
			read_32bit(Block_Memory_Access_Reg_Start_Addr + 4*NofBlocks*BlockSize, last_data_block, NofLastBlock);
			fwrite(last_data_block, sizeof(u32), NofLastBlock, filePointer);
		}
	}
	else if ((NofWords > 8192) && (NofWords <= 16384))
	{
		//First piece of Memory with Block Select Regiter = 0
		write_32bit(Memory_Block_Select, All_Bits_Low);
		for (int i=0; i<81; i++)
		{
			u32 data_word[100];
			read_32bit(Block_Memory_Access_Reg_Start_Addr+4*100*i, data_word, 100);
			fwrite(data_word, sizeof(u32), 100, filePointer);
		}
		u32 data_word[92];
		read_32bit(Block_Memory_Access_Reg_Start_Addr+4*81*100, data_word, 92);
		fwrite(data_word, sizeof(u32), 92, filePointer);

		//Second piece of Memory with Block Select Regiter = 1
		write_32bit(Memory_Block_Select, 1);
		int NofBlocks = ((NofWords - 8192) / BlockSize);
		int NofLastBlock = ((NofWords - 8192) % BlockSize);
		for (int i=0; i<NofBlocks; i++)
		{
			u32 data_word[BlockSize];
			read_32bit(Block_Memory_Access_Reg_Start_Addr+4*BlockSize*i, data_word, BlockSize);
			fwrite(data_word, sizeof(u32), BlockSize, filePointer);
		}
		if(NofLastBlock != 0)
		{
			u32 last_data_block[NofLastBlock];
			read_32bit(Block_Memory_Access_Reg_Start_Addr + 4*NofBlocks*BlockSize, last_data_block, NofLastBlock);
			fwrite(last_data_block, sizeof(u32), NofLastBlock, filePointer);
		}
	}
	else if ((NofWords > 16384) && (NofWords <= 20000))
	{
		//First and Second piece of Memory with Block Select Regiter = 0 and 1
		for (int j=0; j<2; j++)
		{
			write_32bit(Memory_Block_Select, j);
			for (int i=0; i<81; i++)
			{
				u32 data_word[100];
				read_32bit(Block_Memory_Access_Reg_Start_Addr+4*100*i, data_word, 100);
				fwrite(data_word, sizeof(u32), 100, filePointer);
			}
			u32 data_word[92];
			read_32bit(Block_Memory_Access_Reg_Start_Addr+4*81*100, data_word, 92);
			fwrite(data_word, sizeof(u32), 92, filePointer);
		}

		//Third piece of Memory with Block Select Regiter = 2
		write_32bit(Memory_Block_Select, 2);
		int NofBlocks = ((NofWords - 16384) / BlockSize);
		int NofLastBlock = ((NofWords - 16384) % BlockSize);
		for (int i=0; i<NofBlocks; i++)
		{
			u32 data_word[BlockSize];
			read_32bit(Block_Memory_Access_Reg_Start_Addr+4*BlockSize*i, data_word, BlockSize);
			fwrite(data_word, sizeof(u32), BlockSize, filePointer);
		}
		if(NofLastBlock != 0)
		{
			u32 last_data_block[NofLastBlock];
			read_32bit(Block_Memory_Access_Reg_Start_Addr + 4*NofBlocks*BlockSize, last_data_block, NofLastBlock);
			fwrite(last_data_block, sizeof(u32), NofLastBlock, filePointer);
		}
	}
	else 
	{
		printf("No events has been detected.\n");
	}

	fclose(filePointer);

	// Reset the Block Memory
	write_32bit(Readout_Done, All_Bits_High);

    }
	return 0;
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
	for (int i=0; i<64; i++)
	{
		Counters->IO_rate_Counter[i] = tmp_IO[i];
		Counters->IO_total_counts += tmp_IO[i];
	}
}

void Log_Counters(struct Counters_Struct *Counters)
{
	printf("Here are the counters:\n");

	printf("Clock: %llu\n", Counters->ClockCounter);
	printf("Event: %llu\n", Counters->EventCounter);

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

	printf("TB Busy: %llu\n", Counters->TB_BusyCounter);
	printf("CoBo Busy: %llu\n", Counters->CoBo_BusyCounter);
	printf("Transit Busy: %llu\n", Counters->Transit_BusyCounter);
	printf("Buffer Busy: %llu\n", Counters->Buffer_BusyCounter);
	printf("Memory Full: %llu\n", Counters->Mem_FullCounter);
	printf("Memory Reading: %llu\n", Counters->Mem_ReadingCounter);
	printf("Memory Writing: %llu\n", Counters->Mem_WritingCounter);

	for (int i=0; i<64; i++)
	{
		printf("IO#%d: %d\n",Counters->IO_rate_Counter[i]);
	}
	printf("Total IO Counts: %d\n",Counters->IO_total_counts);
}
