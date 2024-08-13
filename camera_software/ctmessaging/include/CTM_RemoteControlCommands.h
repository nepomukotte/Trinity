#ifndef RCTRL_H
#define RCTRL_H

#define NREMOTECTRLCMDS			28
#define TIME_INTERVALS			6

#define REBOOT_CPU_ID			0x01
#define SHUTDOWN_CPU_ID			0x02
#define Init_SIAB_ID			0x03
#define Init_TRGB_ID			0x04
#define Init_COBO_ID			0x05
#define Init_LVPS_ID			0x06
#define Init_PDU_ID				0x07
#define Flush_FTP_ID			0x08

#define UPDATE_RUNTIME		0x09
#define RUN_START				0x0A
#define RUN_STOP				0x0B
#define STATE_MSG_EN_ID			0x0C
#define STATE_MSG_DIS_ID		0x0D
#define GET_LOG_ID				0x0E
#define DO_TRG_SCAN				0x0F

#define Kill_SIAB_ID			0x10
#define Kill_TRGB_ID			0x11
#define Kill_COBO_ID			0x12
#define Kill_LVPS_ID			0x13
#define Kill_PDU_ID				0x14
#define Kill_dataRouter_ID		0x15
#define Kill_getEccServer_ID	0x16
#define Kill_ALL_ID				0x17

#define Fix_LVPS_Client_ID		0x18
#define Archive_LOGS_ID			0x19

#define DO_TRG_SCAN_SF			0x1A
#define DOWNLOAD_STARLINK			0x1B
#define HV_SCAN					0x1C

typedef struct{
	const char *names[NREMOTECTRLCMDS]; 
	unsigned short values[NREMOTECTRLCMDS]; 
} remDict;

typedef struct{ 
	int values[6]; 
} tInterval;

remDict remoteControlCommands = { .names = {"CPUreboot\0",
					"CPUshutdown\0",
					"initSIAB\0",
					"initTriggerBoard\0",
					"initCoBo\0",
					"initLVPS\0",
					"initPDU\0",
					"flushFTP\0",
					"updateRuntime\0",
					"startDAQ\0",
					"stopDAQ\0",
					"SMenable\0",
					"SMdisable\0",
					"getLOGs\0",
					"TriggerScanBF\0",
					"killSIAB\0",
					"killTriggerBoard\0",
					"killCoBo\0",
					"killLVPS\0",
					"killPDU",
					"killdataRouter\0",
					"killECCServer\0",
					"killAll\0",
					"fixLVPS\0",
					"archiveLOGs\0",
					"TriggerScanSF\0",
					"DownloadStarlink\0",
					"HVScan\0"},
				.values = {REBOOT_CPU_ID,
					SHUTDOWN_CPU_ID,
					Init_SIAB_ID,
					Init_TRGB_ID,	
					Init_COBO_ID,	
					Init_LVPS_ID,
					Init_PDU_ID,
					Flush_FTP_ID,
					UPDATE_RUNTIME,
					RUN_START,
					RUN_STOP,
					STATE_MSG_EN_ID,
					STATE_MSG_DIS_ID,
					GET_LOG_ID,
					DO_TRG_SCAN,
					Kill_SIAB_ID,
					Kill_TRGB_ID,
					Kill_COBO_ID,
					Kill_LVPS_ID,
					Kill_PDU_ID,
					Kill_dataRouter_ID,
					Kill_getEccServer_ID,
					Kill_ALL_ID,
					Fix_LVPS_Client_ID,
					Archive_LOGS_ID,
					DO_TRG_SCAN_SF,
					DOWNLOAD_STARLINK,
					HV_SCAN
					}};


tInterval timeInterval = {
	.values = {1,2,5,10,20,30}

};

#endif