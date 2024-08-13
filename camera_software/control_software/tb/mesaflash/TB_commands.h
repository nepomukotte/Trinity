#ifndef __COMMANDS_H
#define __COMMANDS_H

#define CMD_PING			0x01
#define CMD_INIT			0x02
#define CMD_PROGRAM			0x03
#define CMD_PREPARE			0x04
#define CMD_CONFIG			0x05
#define CMD_READ			0x06
#define CMD_WRITE			0x07
#define CMD_START			0x08
#define CMD_STOP			0x09
#define CMD_GL_CTR			0x0A
#define CMD_CL_CTR			0x0B
#define CMD_SAVE			0x0C
#define CMD_CLEAR_COUNTERS	0x0D

static board_access_t access_board;

int read_32bit(u16 address, u32 data[], u8 size);
int write_32bit(u16 address, u32 data);
void En_Trig_Type(u32 trig_type);
void En_BF_Trig(int strech_length_code);
void En_Test_Trig(int strech_length_code, int prescale_value);
void En_Int_Trig(int trig_rate_code, int led_mode);
void En_Ext_Trig(int led_mode);
void En_GPS_Trig(int led_mode);

void Ping_TB(uint32_t* cmd_array, char *response);
void Init_TB(board_access_t *access_board, uint32_t* cmd_array);
int  Program_TB(board_access_t *access_board, uint32_t* cmd_array);
void Clear_All_Counters(uint32_t* cmd_array);
void Prepare_TB(uint32_t* cmd_array);
void Config_TB(uint32_t* cmd_array);
int  Read_Address(uint32_t* cmd_array);
int  Write_Address(uint32_t* cmd_array);
void Start_Trigger(uint32_t* cmd_array);
void Stop_Trigger(uint32_t* cmd_array);
void Get_Global_Ctrs(uint32_t* cmd_array);
void Get_Channel_Ctrs(uint32_t* cmd_array);
void  Save_Events(uint32_t* cmd_array);
board_access_t Process_CMD(uint32_t* cmd_array, char *response, board_access_t ab);

#endif
