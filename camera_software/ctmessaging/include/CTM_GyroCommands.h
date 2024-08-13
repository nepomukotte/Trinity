#ifndef GYRO_H
#define GYRO_H

#define NGYROCOMMANDS 4

#define CMD_G_CONFIG	0x26
#define CMD_G_ENABLE	0x27
#define CMD_G_DISABLE	0x28
#define CMD_G_ARXIV		0x29

typedef struct{
	const char *names[NGYROCOMMANDS]; 
	int values[NGYROCOMMANDS]; 
} gyroDict;

gyroDict gyroCommands = { .names = {"config\0",
					"enable\0",
					"disable\0",
					"save\0"},
				.values = {CMD_G_CONFIG,
					CMD_G_ENABLE,
					CMD_G_DISABLE,
					CMD_G_ARXIV
					}};
#endif