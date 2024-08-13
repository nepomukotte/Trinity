#ifndef _CONSTANTS_H
#define _CONSTANTS_H

#define QID_RC		0x10
#define QID_SIAB	0x20
#define QID_TRGB	0x30
#define QID_COBO	0x40
#define QID_LVPS	0x50
#define QID_PDU		0x60

#define REBOOT_CPU_ID		0x01
#define SHUTDOWN_CPU_ID		0x02
#define Init_SIAB_ID		0x03
#define Init_TRGB_ID		0x04
#define Init_COBO_ID		0x05
#define Init_LVPS_ID		0x06
#define Init_PDU_ID			0x07
#define Flush_FTP_ID		0x08

#define Update_Runtime_ID	0x09
#define RUN_START			0x0A
#define RUN_STOP			0x0B
#define STATE_MSG_EN_ID		0x0C
#define STATE_MSG_DIS_ID	0x0D
#define GET_LOG_ID			0x0E
#define DO_TRG_SCAN			0x0F

#define Kill_SIAB_ID		0x10
#define Kill_TRGB_ID		0x11
#define Kill_COBO_ID		0x12
#define Kill_LVPS_ID		0x13
#define Kill_PDU_ID			0x14
#define Kill_ALL_ID			0x15


#define Process_OFF			0x32
#define Process_ON			0x40

const std::string CS_DIR = "/home/cherenkov/Programs/control_software/";
const std::string RUN_DIR = "/home/cherenkov/Documents/Palestine_Tests/RawDataUnmerged/";
const std::string LOG_DIR = CS_DIR+"fcutils/test/LOGS/";
const std::string BACKUP_DIR = "/home/cherenkov/Backup/";
const std::string FTP_DIR = "/srv/ftp/data/";

// Full Commands for Controlling Trigger Board
const std::string TB_Init					= "0200000000000000";
const std::string TB_Prepare				= "0400000000000000";
const std::string TB_Config_BF_Only			= "05010002000F0003";
const std::string TB_Tscan_Config			= "0501010200005002";	// This one is not using cobo busy signal, in case it was off.
const std::string TB_Config_BF_IntHLED1		= "05050002000F51C3";
const std::string TB_Config_BF_IntHLED10	= "05050002000F54C3";
const std::string TB_Enable_Trigger			= "0800000000000000";
const std::string TB_Disable_Trigger		= "0900000000000000";
const std::string TB_Get_GL_ctrs   			= "0A00000000000000";
const std::string TB_Get_CL_ctrs   			= "0B00000000000000";
const std::string TB_Save_Data				= "0C00000000000000";
const std::string TB_Clear_Counters			= "0D00000000000000";

const std::string CoBo_Start_Run 	 		= "0400000000000000";
const std::string CoBo_Stop_Run 	 		= "0500000000000000";
const std::string CoBo_ReConfigure	 		= "0600000000000000";

const std::string GET_ALL_SIAB_SIPM_TEMP	= "FFFFFFFF0A000000";
const std::string GET_ALL_SIAB_UC_TEMP		= "FFFFFFFF0B000000";
const std::string GET_ALL_SIAB_HV_STATE		= "FFFFFFFF21000000";

#endif