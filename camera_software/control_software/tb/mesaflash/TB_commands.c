#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
#include "anyio.h"
#include "lbp16.h"
#include "eth_boards.h"
#include "TB_Map.h"
#include "TB_commands.h"


u32 trig_rate_list[16] = {Trig_Rate_0Hz,   Trig_Rate_1Hz,  Trig_Rate_3Hz,  Trig_Rate_5Hz,   Trig_Rate_10Hz,
						  Trig_Rate_20Hz,  Trig_Rate_30Hz, Trig_Rate_50Hz, Trig_Rate_100Hz, Trig_Rate_200Hz,
						  Trig_Rate_500Hz, Trig_Rate_1kHz, Trig_Rate_2kHz, Trig_Rate_10kHz, Trig_Rate_100kHz,
						  Trig_Rate_1MHz};

u32 strech_length_list[8] = {strech_length_10ns,  strech_length_20ns,  strech_length_50ns,  strech_length_100ns,
							 strech_length_200ns, strech_length_300ns, strech_length_500ns, strech_length_1us};

u32 disc_deadtime_list[16] = {Deadtime_10ns,  Deadtime_20ns,  Deadtime_30ns,  Deadtime_40ns,  Deadtime_50ns,
							  Deadtime_100ns, Deadtime_150ns, Deadtime_200ns, Deadtime_300ns, Deadtime_400ns,
							  Deadtime_500ns, Deadtime_600ns, Deadtime_700ns, Deadtime_800ns, Deadtime_900ns,
							  Deadtime_1us};

u32 counter_period_list[6] = {One_Second_Interval, Two_Second_Interval, Five_Second_Interval, Ten_Second_Interval, Twenty_Second_Interval, Thirty_Second_Interval};

int read_32bit(u16 address, u32 data[], u8 size)
{
	int send = 0; int recv = 0;
	lbp16_cmd_addr packet;
	u32 command = ((0x4200) | (0x80 | size));
	LBP16_INIT_PACKET4(packet, command, address);
	send = eth_socket_send_packet(&packet, sizeof(packet));
	usleep(25*1000);
	recv = eth_socket_recv_packet(data, size*sizeof(data[0]));
	return recv;
}

int write_32bit(u16 address, u32 data)
{
	int send = 0;
	u32 tmp[1] = {0};
	lbp16_cmd_addr_data32 packet;
	LBP16_INIT_PACKET8(packet, 0xC281, address, data);
	printf("Adress to write: %#04x\n", address);
	printf("Data supposed to write: %lu\n", data);
	send = eth_socket_send_packet(&packet, sizeof(packet));
	return send;
}

void En_BF_Trig(int strech_length_code)
{
	write_32bit(Disc_Stretch_Enable_0, All_Bits_High);
	write_32bit(Disc_Stretch_Enable_1, All_Bits_High);
	write_32bit(Disc_Stretch_Length, strech_length_list[strech_length_code]);
}

void En_Test_Trig(int strech_length_code, int prescale_value)
{
	write_32bit(Disc_Stretch_Enable_0, All_Bits_High);
	write_32bit(Disc_Stretch_Enable_1, All_Bits_High);
	write_32bit(Disc_Stretch_Length, strech_length_list[strech_length_code]);
	write_32bit(Disc_Test_Trigger_Prescale, prescale_value);
}

void En_Int_Trig(int trig_rate_code, int led_mode)
{
	write_32bit(Internal_Trigger_Prescale, trig_rate_list[trig_rate_code]);

	if(led_mode == 0) {
		write_32bit(Internal_Trigger_Mode, All_Bits_Low);
	}else if(led_mode == 1){
		write_32bit(Internal_Trigger_Mode, First_Bit_True);
	}else{
		printf("Could not understand the LED Mode option for internal trigger\n");
	}
}

void En_Ext_Trig(int led_mode)
{
	if(led_mode == 0) {
		write_32bit(External_Trigger_Mode, All_Bits_Low);
	}else if(led_mode == 1){
		write_32bit(External_Trigger_Mode, First_Bit_True);
	}else{
		printf("Could not understand the LED Mode option for external trigger\n");
	}
}

void En_GPS_Trig(int led_mode)
{
	if(led_mode == 0) {
		write_32bit(GPS_Trigger_Mode, All_Bits_Low);
	}else if(led_mode == 1){
		write_32bit(GPS_Trigger_Mode, First_Bit_True);
	}else{
		printf("Could not understand the LED Mode option for GPS trigger\n");
	}
}

void Ping_TB(uint32_t* cmd_array, char *response)
{
	strcpy(response, "3333");
}

void Init_TB(board_access_t *access_board, uint32_t* cmd_array)
{
	printf("---------- Starting to Initialize TB ------------\n");
	access_board->eth = 1;
	access_board->address = 1;
	access_board->dev_addr = "192.168.3.121";
	access_board->device_name = "7I80";


	if(anyio_list_dev(access_board)==0){
		eth_boards_init(access_board);
		eth_boards_scan(access_board);
	}
	else{
		printf("Trigger Board has been already initialized\n");
	}
	//anyio_list_dev(access_board);
	printf("---------- Finished Initializing TB -----------\n\n");
}

int Program_TB(board_access_t *access_board, uint32_t* cmd_array)
{
	printf("---------- Starting to Upload TB Logic --------\n");
	int ret = -1;
	board_t *board = NULL;
	static int fallback_flag = 0;
	static char bitfile_name[255];
	char bitfile_str[] = "April_1st_Version.bit";
	strncpy(bitfile_name, bitfile_str, sizeof(bitfile_name));
	board = anyio_get_dev(access_board);
	if(board == NULL)
	{
		printf("No board was found\n");
		return -1;
	}
	board->open(board);
	ret = anyio_dev_write_flash(board, bitfile_name, fallback_flag);
	board->close(board);
	anyio_cleanup(access_board);
	printf("---------- Finished Uploading TB Logic --------\n\n");
	return ret;
}

void Prepare_TB(uint32_t* cmd_array)
{
	printf("---------- Starting to Prepare TB -------------\n");

	uint32_t Transit_Busy_Duration = (((cmd_array[1] << 8) | cmd_array[2]) & 0xFFFF);

	// Clearing the Busy Bits
	printf("  Clearing the Busy Registers\n");
	write_32bit(Clear_Busy, First_Bit_True);

	// Enable the Block Memory
	printf("\n  Enabling the Block Memory\n");
	write_32bit(Block_Memory_Enable, First_Bit_True);

	// Set the Block Memory Start Address to zero
	printf("\n  Setting the Block Memory Address to zero\n");
	write_32bit(Memory_Start_Address, All_Bits_Low);

	// Reset the Block Memory
	printf("\n  Resetting the Block Memory\n");
	write_32bit(Readout_Done, First_Bit_True);

	// Set the Transit Busy Length to (Transit_Busy_Duration*10) ns
	printf("\n  Setting the Transit Busy length to %dns\n", Transit_Busy_Duration*10);
	write_32bit(Transit_Busy_Length, Transit_Busy_Duration);

	printf("---------- Finished Preparing TB --------------\n\n");
}

void Config_TB(uint32_t* cmd_array)
{
	printf("---------- Starting to Configure TB -----------\n");

	int trigger_type = (cmd_array[1] & 0x1F);
	int Bifocal_En   = ((trigger_type & 0x01)>0);
	int DiscTest_En  = ((trigger_type & 0x02)>0);
	int Internal_En  = ((trigger_type & 0x04)>0);
	int External_En  = ((trigger_type & 0x08)>0);
	int GPS_En = ((trigger_type & 0x10)>0);

	int ctr_integ_code = (cmd_array[2] & 0xFF);
	int Disc_Stretch_Length_code = (cmd_array[3] & 0xFF);
	int Disc_Test_Prescale_Value = (cmd_array[4] & 0xFF);
	int HLED_Delay   = (cmd_array[5] & 0xFF);
	int Disc_Dt_code = (cmd_array[6] & 0xF0);

	int Int_Trig_rate_code = (cmd_array[6] & 0x0F);
	int Int_Trig_led_mode  = ((cmd_array[7] & 0xC0)>0);
	int Ext_Trig_led_mode  = ((cmd_array[7] & 0x30)>0);
	int GPS_Trig_led_mode  = ((cmd_array[7] & 0x0C)>0);
	int Enable_Busy_mode   = (cmd_array[7] & 0x03);

	printf("Trigger type: %d\n", trigger_type);
	printf("Bifocal Trigger Flag: %d\n", Bifocal_En);
	printf("Disc Test Trigger Flag: %d\n", DiscTest_En);
	printf("Internal Trigger Flag: %d\n", Internal_En);
	printf("External Trigger Flag: %d\n", External_En);
	printf("GPS Trigger Flag: %d\n", GPS_En);
	printf("Channel Counter Integration Period Code: %d\n", ctr_integ_code);
	printf("Disc Test Trigger Prescale Value: %d\n", Disc_Test_Prescale_Value);
	printf("Discriminator Strech Length Code: %d\n", Disc_Stretch_Length_code);
	printf("HLED_Delay: %d ns\n", HLED_Delay*10);
	printf("Disc_Dt_code: %d\n", Disc_Dt_code);
	printf("Int_Trig_rate_code: %d\n", Int_Trig_rate_code);
	printf("Int_Trig_led_mode: %d\n", Int_Trig_led_mode);
	printf("Ext_Trig_led_mode: %d\n", Ext_Trig_led_mode);
	printf("GPS_Trig_led_mode: %d\n", GPS_Trig_led_mode);
	printf("Enable Busy mode: %d\n", Enable_Busy_mode);

	write_32bit(Rate_Counter_Period, counter_period_list[ctr_integ_code]);
	write_32bit(Clear_Counters, First_Bit_True);
	write_32bit(Disc_Deadtime, disc_deadtime_list[Disc_Dt_code]);

	if(Bifocal_En == 1){
		En_BF_Trig(strech_length_list[Disc_Stretch_Length_code]);
	}

	if(DiscTest_En == 1){
		En_Test_Trig(strech_length_list[Disc_Stretch_Length_code], Disc_Test_Prescale_Value);
	}else{
		write_32bit(Disc_Test_Trigger_Prescale, 0);	// This is to make sure test trigger rate is zero, if it is not being enabled.
	}

	if(Internal_En == 1){
		En_Int_Trig(Int_Trig_rate_code, Int_Trig_led_mode);
	}else{
		En_Int_Trig(0, 0);	// This is to make sure HLED is not flashing from previous configurations.
	}

	if(External_En == 1){
		En_Ext_Trig(Ext_Trig_led_mode);
	}else{
		En_Ext_Trig(0);	// This is to make sure HLED is not flashing from previous configurations.
	}

	if(GPS_En == 1){
		En_GPS_Trig(Ext_Trig_led_mode);
	}else{
		En_GPS_Trig(0);	// This is to make sure HLED is not flashing from previous configurations.
	}

	write_32bit(Enable_Busy, Enable_Busy_mode);
	write_32bit(LED_Delay, HLED_Delay);
	write_32bit(Enable_Trigger_Types, trigger_type);

	printf("---------- Finished Configuring TB ------------\n\n");
}

int Read_Address(uint32_t* cmd_array)
{
	printf("---------- Reading a Register -----------------\n");
	int ret = 0;
	u16 address = (((cmd_array[1] << 8) | cmd_array[2]) & 0xFFFF);
	u8 size = (cmd_array[3] & 0xFF);
	u32 data[size];
	printf("Adress to read: %#04x\n", address);
	printf("Number of words to read: %d\n", size);
	ret = read_32bit(address, data, size);
	for(int i=0; i<size; i++){
		printf("Word#%d: %d\n", i, data[i]);
	}
	printf("---------- Finished Reading a Register --------\n\n");
	return ret;
}

int Write_Address(uint32_t* cmd_array)
{
	printf("---------- Writing to a Register --------------\n");
	int ret = 0;
	u16 address = (((cmd_array[1] << 8) | cmd_array[2]) & 0xFFFF);
	u32 data =    (((cmd_array[3] << 24) | (cmd_array[4] << 16) | (cmd_array[5] << 8) | (cmd_array[4])) & 0xFFFFFFFF);
	printf("Adress to write: %#04x\n", address);
	printf("Value to write: %d\n", data);
	ret = write_32bit(address, data);
	printf("---------- Finished Writing to a Register -----\n\n");
	return ret;
}

void Start_Trigger(uint32_t* cmd_array)
{
	printf("---------- Starting the Trigger ---------------\n");
	write_32bit(Enable_Triggering, First_Bit_True);
	printf("-----------------------------------------------\n\n");
}

void Stop_Trigger(uint32_t* cmd_array)
{
	printf("---------- Stopping the Trigger ---------------\n");
	write_32bit(Enable_Triggering, All_Bits_Low);
	printf("-----------------------------------------------\n\n");
}

void Get_Global_Ctrs(uint32_t* cmd_array)
{
	printf("---------- Starting to Record the GC ----------\n");
	// Reading 29 registers starting from Counters_First_Addr (0x2020)
	// This array will be taken apart to calculate individual counters
	u32 tmp[29];
	u64 counters[21];
	u32 counters_overflow[3];
	
	u8 counter_type = (cmd_array[1] & 0xFF);
	u32 counter_run_num = (((cmd_array[2] << 8) | cmd_array[3]) & 0xFFFF);

	write_32bit(Save_Counters, First_Bit_True);
	read_32bit(Counters_First_Addr, tmp, 29);
	read_32bit(Overflow_0, counters_overflow, 3);

	counters[0]  = (((u64)tmp[1] << 32 )| (u64)tmp[0]);
	counters[1]  = tmp[2];
	counters[2]  = tmp[3];
	counters[3]  = tmp[4];
	counters[4]  = tmp[5];
	counters[5]  = tmp[6];
	counters[6]  = tmp[7];
	counters[7]	 = tmp[8];
	counters[8]	 = tmp[9];
	counters[9]	 = tmp[10];
	counters[10] = tmp[11];
	counters[11] = tmp[12];
	counters[12] = tmp[13];
	counters[13] = tmp[14];
	counters[14] = (((u64)tmp[16] << 32 )| (u64)tmp[15]);
	counters[15] = (((u64)tmp[18] << 32 )| (u64)tmp[17]);
	counters[16] = (((u64)tmp[20] << 32 )| (u64)tmp[19]);
	counters[17] = (((u64)tmp[22] << 32 )| (u64)tmp[21]);
	counters[18] = (((u64)tmp[24] << 32 )| (u64)tmp[23]);
	counters[19] = (((u64)tmp[26] << 32 )| (u64)tmp[25]);
	counters[20] = (((u64)tmp[28] << 32 )| (u64)tmp[27]);

	printf("\nHere are the counters:\n");
	printf("Clock: %lu\n", counters[0]);
	printf("Event: %u\n", counters[1]);
	printf("Bifocal: %d\n", counters[2]);
	printf("Disc Test: %d\n", counters[3]);
	printf("Internal: %d\n", counters[4]);
	printf("External: %d\n", counters[5]);
	printf("GPS: %d\n", counters[6]);
	printf("HLED: %d\n", counters[7]);
	
	printf("Bifocal Live: %d\n", counters[8]);
	printf("Disc Test Live: %d\n", counters[9]);
	printf("Internal Live: %d\n", counters[10]);
	printf("External Live: %d\n", counters[11]);
	printf("GPS Live: %d\n", counters[12]);
	printf("HLED Live: %d\n", counters[13]);

	printf("TB Busy: %lu\n", counters[14]);
	printf("CoBo Busy: %lu\n", counters[15]);
	printf("Transit Busy: %lu\n", counters[16]);
	printf("Buffer Busy: %lu\n", counters[17]);
	printf("Memory Full: %lu\n", counters[18]);
	printf("Memory Reading: %lu\n", counters[19]);
	printf("Memory Writing: %lu\n", counters[20]);

	FILE *filePointer;
	char filename[200];
	sprintf(filename, "/home/cherenkov/Programs/control_software/tb/mesaflash/G_Ctrs_%01d_%05d.bin", counter_type, counter_run_num);
	filePointer = fopen(filename, "a");
	fwrite(counters, sizeof(u64), 21, filePointer);
	fclose(filePointer);
	printf("Global Counters Filename: %s\n", filename);
	printf("---------- Finished Recording the GC ----------\n\n");
}

void Get_Channel_Ctrs(uint32_t* cmd_array)
{
	printf("---------- Starting to Record the CC ----------\n");
	u32 tmp_IO[64];
	u8 counter_type = (cmd_array[1] & 0xFF);
	u32 counter_run_num = (((cmd_array[2] << 8) | cmd_array[3]) & 0xFFFF);

	write_32bit(Save_Counters, First_Bit_True);
	read_32bit(IO_0_rate_Addr, tmp_IO, 64);

	printf("\n");
	for (int i=0; i<64; i++)
	{
		printf("IO#%d: %d\n", i, tmp_IO[i]);
	}

	FILE *filePointer;
	char filename[100];
	sprintf(filename, "/home/cherenkov/Programs/control_software/tb/mesaflash/C_Ctrs_%01d_%05d.bin", counter_type, counter_run_num);
	filePointer = fopen(filename, "a");
	fwrite(tmp_IO, sizeof(u32), 64, filePointer);
	fclose(filePointer);
	printf("Channel Counters Filename: %s\n", filename);
	printf("---------- Finished Recording the CC ----------\n\n");
}

void Save_Events(uint32_t* cmd_array)
{
	printf("---------- Saving Events to the Disk ----------\n");
	FILE *filePointer;
	u32 NofWords = 0;
	int BlockSize = 100;
	u32 NofEvents[1] = {0};
	int Run_Number = (((cmd_array[1] << 8) | cmd_array[2]) & 0xFFFF);

	// Starting the Readout Process
	write_32bit(Readout_Start, First_Bit_True);
	sleep_ns(5*1000*1000);

	// Reading the number of events in the Memory until it makes sense
	int tol_it = 10;
	int it_ctr = 0;
	do{
		read_32bit(Events_Written, NofEvents, 1);
		printf("No. Events at Loop %d: %d\n",it_ctr,NofEvents[0]);
		it_ctr++;
	}while(NofEvents[0]>4000 && it_ctr<tol_it);
	printf("Attempts To Read No. of Events: %d\n",it_ctr);

	NofWords = (5 * NofEvents[0]);
	printf("Number of Words to read: %d\n", NofWords);

	// Writing Events_Writtents to a file
	char filename[200];
	char run_num_tail[10];
	struct tm *timenow;
	time_t now = time(NULL);
	timenow = localtime(&now);
	strftime(filename, sizeof(filename), "/home/cherenkov/Programs/control_software/tb/mesaflash/TB_data_%Y_%m_%d_%H_%M_%S_", timenow);
	sprintf(run_num_tail, "%04d.bin", Run_Number);
	strcat(filename, run_num_tail);
	filePointer = fopen(filename, "wb");

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
	write_32bit(Readout_Done, First_Bit_True);

	printf("---------- Finished Saving the Events ---------\n\n");
}

void Clear_All_Counters(uint32_t* cmd_array)
{
	printf("---------- Starting to Clear TB Counters ------\n");
	write_32bit(Clear_Counters, First_Bit_True);
	printf("---------- Finished Clearing TB Counters ------\n\n");
}

board_access_t Process_CMD(uint32_t* cmd_array, char *response, board_access_t ab)
{
	uint32_t CMD_ID = cmd_array[0];
	printf("Command ID is: %d\n", CMD_ID);

	access_board = ab;

	switch (CMD_ID)
	{
		case CMD_PING:
			Ping_TB(cmd_array, response);
			break;
		case CMD_INIT:
			Init_TB(&access_board, cmd_array);
			break;
		case CMD_PROGRAM:
			Program_TB(&access_board, cmd_array);
			break;
		case CMD_PREPARE:
			Prepare_TB(cmd_array);
			break;
		case CMD_CONFIG:
			Config_TB(cmd_array);
			break;
		case CMD_READ:
			Read_Address(cmd_array);
			break;
		case CMD_WRITE:
			Write_Address(cmd_array);
			break;
		case CMD_START:
			Start_Trigger(cmd_array);
			break;
		case CMD_STOP:
			Stop_Trigger(cmd_array);
			break;
		case CMD_GL_CTR:
			Get_Global_Ctrs(cmd_array);
			break;
		case CMD_CL_CTR:
			Get_Channel_Ctrs(cmd_array);
			break;
		case CMD_SAVE:
			Save_Events(cmd_array);
			break;
		case CMD_CLEAR_COUNTERS:
			Clear_All_Counters(cmd_array);
			break;
		default:
			break;
	}

	return access_board;
}
