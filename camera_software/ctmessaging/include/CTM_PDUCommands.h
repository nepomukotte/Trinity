#ifndef PDU_H
#define PDU_H

#define NPDUCOMMANDS 6
#define NPDUDEVICES 2

#define CMD_SETUP		0x01
#define CMD_TURN_ON		0x02
#define CMD_TURN_OFF 	0x03
#define CMD_GET_STATE	0x04
#define CMD_GET_USAGE	0x05

#define CMD_PING_PDU	0x30

typedef struct{
	const char *names[NPDUCOMMANDS]; 
	int values[NPDUCOMMANDS]; 
} pduDict;

typedef struct{
	const char *names[NPDUDEVICES];
	int values[NPDUDEVICES];
} pduDev;

pduDict pduCommands = { .names = {"init\0",
					"powerOn\0",
					"powerOff\0",
					"getState\0"
					"getUsage\0"
					"ping\0"},
				.values = {CMD_SETUP,
					CMD_TURN_ON,
					CMD_TURN_OFF,
					CMD_GET_STATE,
					CMD_GET_USAGE,
					CMD_PING_PDU	
					}};

pduDev pduDevices = {.names={"lvps","pump"},
					.values={9,8}};
#endif