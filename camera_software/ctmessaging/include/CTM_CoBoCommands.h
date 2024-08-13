#ifndef COBO_H
#define COBO_H

#define NCOBOCOMMANDS 7

#define START_ECC_SERVER	0x01
#define START_DATA_ROUTER	0x02
#define INIT_EXPERIMENT		0x03
#define START_RUN		0x04
#define STOP_RUN		0x05
#define AsAd_OFF		0x06

#define CMD_PING_COBO		0x30

typedef struct{
	const char *names[NCOBOCOMMANDS]; 
	const int values[NCOBOCOMMANDS];
} coboDict;


coboDict coboCommands = { .names = {"startEccServer\0",
					"startDataRouter\0",
					"initExperiment\0",
					"startDAQ\0",
					"stopDAQ\0",
					"stopASAD\0",
					"ping\0"},
				.values = {START_ECC_SERVER,
					START_DATA_ROUTER,
					INIT_EXPERIMENT,
					START_RUN,
					STOP_RUN,
					AsAd_OFF,
					CMD_PING_COBO
					}};
#endif