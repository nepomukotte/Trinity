#ifndef PROC_H
#define PROC_H

#define NTOTALPROCESS			8
#define QID_RC					0x10
#define QID_SIAB				0x20
#define QID_GYRO				0x21
#define QID_TRGB				0x30
#define QID_COBO				0x40
#define QID_LVPS				0x50
#define QID_PDU					0x60
#define QID_SEQ					0xAA

typedef struct{
	const char* names[NTOTALPROCESS]; 
	int values[NTOTALPROCESS]; 
} procIDDict;

procIDDict processID = {.names = {"RemoteControl\0",
					"SIAB\0",
					"TriggerBoard\0",
					"CoBo\0",
					"LVPS\0",
					"PDU\0",
					"Sequence\0",
					"Gyro\0"}, 
					.values = {QID_RC,
					QID_SIAB,
					QID_TRGB,
					QID_COBO,
					QID_LVPS,
					QID_PDU,
					QID_SEQ,
					QID_GYRO}};
#endif