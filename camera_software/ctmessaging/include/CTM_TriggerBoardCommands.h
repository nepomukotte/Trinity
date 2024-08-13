#ifndef TRGB_H
#define TRGB_H

#define NTBCOMMANDS			14
#define NTTYPES 5
#define NDISCCODE 8
#define NDEADCODE 16
#define NTRIGRATECODE 16
#define NBUSYCODE 4

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
#define CMD_SAVE_COUNTERS	0x0E

typedef struct{
	const char *names[14]; 
	int values[14]; 
} tbDict;

typedef struct{
	const char *names[5]; 
	int values[5]; 
} tbTrigTypes;

typedef struct{
	int values[8]; 
} tbDiscCode;

typedef struct{
	int values[16]; 
} tbDeadTimeCode;

typedef struct{
	int values[16]; 
} tbIntTrigCode;
typedef struct{
	const char *names[4]; 
	int values[4]; 
} tbBusySignal;

tbDict triggerBoardCommands = { .names = {"ping\0",
					"init\0",
					"program\0",
					"prepare\0",
					"config\0",
					"read\0",
					"write\0",
					"start\0",
					"stop\0",
					"globalCounters\0",
					"channelCounters\0",
					"save\0",
					"clearCounters\0",
					"saveCounters\0"},
				.values = {CMD_PING,			
					CMD_INIT,	
					CMD_PROGRAM,	
					CMD_PREPARE,
					CMD_CONFIG,
					CMD_READ,
					CMD_WRITE,
					CMD_START,
					CMD_STOP,
					CMD_GL_CTR,
					CMD_CL_CTR,
					CMD_SAVE,
					CMD_CLEAR_COUNTERS,
					CMD_SAVE_COUNTERS,
					}};

tbTrigTypes triggerBoardTT = {.names = {"bifocal\0",
					"discTest\0",
					"internal\0",
					"external\0",
					"gps\0"},
				.values = {0x01,
				0x02,
				0x04,
				0x08,
				0x10}};

tbDiscCode triggerBoardDC = {.values = {10,
						20,
						50,
						100,
						200,
						300,
						500,
						1000}};

tbDeadTimeCode triggerBoardDTC = {.values = {10,20,30,40,50,
									100,200,300,400,500,
									600,700,800,900,1000}};

tbIntTrigCode triggerBoardIC = {.values = {0,1,3,5,10,20,30,50,
									100,200,500,1000,2000,
									10000,100000,1000000}};

tbBusySignal triggerBoardBusy = {.names={"off\0","cobo\0",
								"buffer\0",
								"both\0"},
								.values = {0,1,2,3}};
#endif