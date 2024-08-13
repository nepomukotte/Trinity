#ifndef SEQUENCES_H
#define SEQUENCES_H

#define NSEQUENCECOMMANDS 6
#define NCONFIGSEQUENCES 6

#define POWER_ON_SEQ	0x01
#define INIT_SEQ	0x02
#define CONFIGURE_SEQ	0x03
#define START_SEQ	0x04
#define STOP_SEQ	0x05
#define POWER_OFF_SEQ	0x06


typedef struct{
	const char *names[NSEQUENCECOMMANDS]; 
	int values[NSEQUENCECOMMANDS]; 
} sequenceDict;

typedef struct{
	const char *names[NCONFIGSEQUENCES];
	int values[NCONFIGSEQUENCES];
} configDict;

sequenceDict seqCommands = { .names = {"powerOn\0",
					"init\0",
					"config\0",
					"start\0",
					"stop\0",
					"powerOff\0"},
				.values = {POWER_ON_SEQ,
					INIT_SEQ,
					CONFIGURE_SEQ,
					START_SEQ,
					STOP_SEQ,
					POWER_OFF_SEQ	
					}};

configDict configSequences = { .names = {"bifocalOnly\0",
					"hled1Hz\0",
					"hled10Hz\0",
					"bifocalHled1Hz\0",
					"bifocalHled10Hz\0",
					"singleFocusOnly"},
					.values = {1,2,3,4,5,6}};
#endif