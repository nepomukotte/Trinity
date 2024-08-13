#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <thread>
#include <functional>
#include <sys/reboot.h>
#include <iomanip>
#include "posix_message_queue.hpp"
#include "fcutils/utility/utility.hpp"
#include "fcutils/string/string.hpp"
#include "state_message_cmds.h"

using namespace std;
using fcutils::system::posix_message_queue;

int Run_Number = 0;
int TScanNo = 0;
int run_interval = 60;			// In seconds
int state_msg_interval = 30;	// In seconds	
bool RUN_ENABLED = false;
bool EN_STATE_MSG = false;
bool IsInsideRunTimer = false;
bool DoNotContinueRun = false;
std::string state_response = "1010";

int SiPM_Temp_List[16] = {0,2,4,6,9,11,13,15,16,18,20,22,25,27,29,31};
int UC_Temp_List[16] = {1,3,5,7,8,10,12,14,17,19,21,23,24,26,28,30};
std::string LogFileID[6] = {"11111", "22222", "33333", "44444", "55555", "66666"};
std::string LogFileList[6] = {"CS", "Client_SIAB", "Client_TB", "Client_CoBo", "Client_LVPS", "Client_PDU"};


std::string mqname[8] = {"/CStoRC","/RCtoCS","/writetoRC","/writetoSIAB","/writetoTB","/writetoCoBo","/writetoPS","/writetoPDU"};
posix_message_queue rqcs(mqname[0].c_str());	// from Control Software to Master Control
posix_message_queue wqcs(mqname[1].c_str());	// from Master Control to Control Software
posix_message_queue rqrc(mqname[2].c_str());	// from clients to Master Control
posix_message_queue wqsiab(mqname[3].c_str());	// from Master Control to SIABs
posix_message_queue wqtrgb(mqname[4].c_str());	// from Master Control to TB
posix_message_queue wqcobo(mqname[5].c_str());	// from Master Control to CoBo
posix_message_queue wqlvps(mqname[6].c_str());	// from Master Control to LVPS
posix_message_queue wqpdu(mqname[7].c_str());	// from Master Control to PDU

// Equivalent Hex of 50 (DEC) is 32 (Hex) OR 2 (Char) showing in state messages. This means process is turned on.
// Equivalent Hex of 64 (DEC) is 40 (Hex) OR @ (Char) showing in state messages. This means process is turned off.
int SIAB_STATE = Process_OFF;
int TRGB_STATE = Process_OFF;
int COBO_STATE = Process_OFF;
int LVPS_STATE = Process_OFF;
int PDU_STATE  = Process_OFF;

std::string Get_DateTime_Str()
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[20];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%y%m%d_%H%M%S",timeinfo);
  std::string str(buffer);
  return str;
}

std::string Get_Date_Str()
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[20];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%y%m%d",timeinfo);
  std::string str(buffer);
  return str;
}

/////////////////////////////////////////////////////////
///////////////// Power-On Seq Commands /////////////////
/////////////////////////////////////////////////////////

void Reset_CT_CPU()
{
	system("sudo reboot &");
}

void Shutdown_CT_CPU()
{
	system("sudo shutdown -h now &");
}

void Init_SIAB()
{
	std::string command = "sudo ./client_siab >> "+LOG_DIR+"Client_SIAB.log &";
	system(command.c_str());
	SIAB_STATE = Process_ON;
}

void Kill_SIAB()
{
	system("sudo echo $(/bin/ps -fu $USER | awk '/client_siab/ && !/awk/ {print $2}' | xargs kill)");
	SIAB_STATE = Process_OFF;
}

void Init_TRGB()
{
	std::string command = "sudo "+CS_DIR+"tb/mesaflash/client_tb >> "+LOG_DIR+"Client_TB.log &";
	system(command.c_str());
	TRGB_STATE = Process_ON;
}

void Kill_TRGB()
{
	system("sudo echo $(/bin/ps -fu $USER | awk '/client_tb/ && !/awk/ {print $2}' | xargs kill)");
	TRGB_STATE = Process_OFF;
}

void Init_COBO()
{
	std::string command = "sudo ./client_cobo >> "+LOG_DIR+"Client_CoBo.log &";
	system(command.c_str());
	COBO_STATE = Process_ON;
}

void Kill_COBO()
{
	system("sudo echo $(/bin/ps -fu $USER | awk '/client_cobo/ && !/awk/ {print $2}' | xargs kill)");
	COBO_STATE =Process_OFF;
}

void Init_LVPS()
{
	std::string command = "sudo ./client_lvps -p /dev/ttyS1 >> "+LOG_DIR+"Client_LVPS.log &";
	system(command.c_str());
	LVPS_STATE = Process_ON;
}

void Kill_LVPS()
{
	system("sudo echo $(/bin/ps -fu $USER | awk '/client_lvps/ && !/awk/ {print $2}' | xargs kill)");
	LVPS_STATE = Process_OFF;
}

void Init_PDU()
{
	std::string command = "sudo ./client_pdu >> "+LOG_DIR+"Client_PDU.log &";
	system(command.c_str());
	PDU_STATE = Process_ON;
}

void Kill_PDU()
{
	system("sudo echo $(/bin/ps -fu $USER | awk '/client_pdu/ && !/awk/ {print $2}' | xargs kill)");
	PDU_STATE = Process_OFF;
}

void Kill_All()
{
	Kill_COBO();
	Kill_PDU();
	Kill_SIAB();
	Kill_TRGB();
	Kill_LVPS();
}

void Flush_FTP()
{
	// Create a backup directory for today's files inside /home/cherenkov/Backup
	std::string Today_DIR = BACKUP_DIR+"Backup_"+Get_Date_Str()+"/";
	std::string Create_CMD = "sudo mkdir "+Today_DIR+" &";
	system(Create_CMD.c_str());
	sleep(1);

	// Make a log of files left in the download folder
	// Place one copy in today's backup directory and one copy in ftp download directory
	std::string Ex_Files_Log = "CT_0_"+Get_DateTime_Str()+"_99999.log";
	std::string List_CMD = "sudo ls -l /srv/ftp/data/ > "+Today_DIR+Ex_Files_Log+ " &";
	system(List_CMD.c_str());
	sleep(1);

	// Move all files left in the download folder to today's backup directory
	std::string Move_CMD = "sudo mv /srv/ftp/data/* "+Today_DIR+" &";
	system(Move_CMD.c_str());
	sleep(60);

	// Place one copy of log of files in ftp download directory
	std::string CopyList_CMD = "sudo cp "+Today_DIR+Ex_Files_Log+" "+FTP_DIR+" &";
	system(CopyList_CMD.c_str());
	sleep(1);

	std::cout << "Finished moving extra files to: " << Today_DIR << std::endl;
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
///////////////// Run Control Commands /////////////////
////////////////////////////////////////////////////////

void Update_Runtime(std::string& msg)
{
	std::string runtime_str = msg.substr(0, 4);
	int runtime = stoi(runtime_str,0,16);
	run_interval = runtime;
}

std::string Update_TB_GC_CMD(int Run_Num)
{
	std::string temp = fcutils::string::int_to_hex(Run_Num);
	std::string run_num_str = temp.substr(6, 4);
	std::string GC_CMD = "0A00" + run_num_str + "00000000";
	return GC_CMD;
}

std::string Update_TB_CC_CMD(int Run_Num)
{
	std::string temp = fcutils::string::int_to_hex(Run_Num);
	std::string run_num_str = temp.substr(6, 4);
	std::string CC_CMD = "0B00" + run_num_str + "00000000";
	return CC_CMD;
}

std::string Update_TB_Event_CMD(int Run_Num)
{
	std::string temp = fcutils::string::int_to_hex(Run_Num);
	std::string run_num_str = temp.substr(6, 4);
	std::string Event_CMD = "0C" + run_num_str + "0000000000";
	return Event_CMD;
}

void Start_Run()
{
	std::cout << "\n---------------------" << std::endl;
	std::cout << "Starting the run#" << Run_Number << std::endl;
	std::cout << "Please wait 5 seconds before stopping the run" << std::endl;
	wqcobo.send(CoBo_Start_Run);
	sleep(5);
	wqtrgb.send(TB_Enable_Trigger);
	sleep(5);
	RUN_ENABLED = true;
	std::cout << "You can stop the run now." << std::endl;
}

void Stop_Triggering()
{
	std::cout << "Stopping the run#" << Run_Number << std::endl;

	std::cout << "Disabling the Global Trigger" << std::endl;
	wqtrgb.send(TB_Disable_Trigger);
	usleep(500*1000);
	std::cout << "Stopping the CoBo ..." << std::endl;
	wqcobo.send(CoBo_Stop_Run);
	usleep(500*1000);
	std::cout << "Saving the Global Counters ..." << std::endl;
	wqtrgb.send(Update_TB_GC_CMD(Run_Number));
	sleep(1);
	//std::cout << "Saving the Channel Counters ..." << std::endl;
	//wqtrgb.send(Update_TB_CC_CMD(Run_Number));
	//sleep(1);
	std::cout << "Saving the Events data ..." << std::endl;
	wqtrgb.send(Update_TB_Event_CMD(Run_Number));
	sleep(2);
	std::cout << "Clearing the Counters ..." << std::endl;
	wqtrgb.send(TB_Clear_Counters);

	Run_Number++;

	std::cout << "Moving the run files to this directory: " << RUN_DIR << std::endl;
	std::string Move_Files = "/home/cherenkov/Programs/control_software/fcutils/test/scripts/MoveFilesDone.sh " + RUN_DIR + " &";
	system(Move_Files.c_str());
	std::cout << "---------------------" << std::endl;
}

void Stop_Run()
{
	if(RUN_ENABLED)
	{
		while(IsInsideRunTimer){
			std::cout << "Waiting 2 Seconds for Run Timer to finish ..." << std::endl;
			sleep(2);
		}
		std::cout << "Disabling Run_Enabled Flag ..." << std::endl;
		RUN_ENABLED = false;
		std::cout << "Stopping the Run Manually ..." << std::endl;
		sleep(2);
		Stop_Triggering();
		sleep(1);
		std::cout << "Run is stopped" << std::endl;
		std::cout << "---------------------" << std::endl;
	}
}

void Check_Run_Status()
{
	if(RUN_ENABLED)
	{
		IsInsideRunTimer = true;
		std::cout << "Generating Run File" << std::endl;
		Stop_Triggering();
		sleep(1);
		Start_Run();			
		IsInsideRunTimer = false;
	}
}

void Get_LOG(std::string& msg)
{
	int Client_ID = stoi(msg.substr(0, 2), 0 ,16);
	int NofLines = stoi(msg.substr(2, 4), 0, 16);

	std::string Out_filename = "CT_0_"+Get_DateTime_Str()+"_"+LogFileID[Client_ID]+".log";
	std::string command = "tail -"+std::to_string(NofLines)+" "+LOG_DIR+LogFileList[Client_ID]+".log > "+Out_filename+" &";
	system(command.c_str());

	std::string Move_File = "sudo mv "+Out_filename+" /srv/ftp/data/ &";
	system(Move_File.c_str());
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
///////////////// State MSG Commands ///////////////////
////////////////////////////////////////////////////////

void Enable_State_MSG()
{
	EN_STATE_MSG = true;
}

void Disable_State_MSG()
{
	EN_STATE_MSG = false;
}

std::string Get_Camera_HV_State()
{
    std::cout << "Getting HV State from SIABs" << std::endl;
    wqsiab.send(GET_ALL_SIAB_HV_STATE);
    sleep(2);

    FILE *hFile = NULL;
    uint32_t *hBuffer = NULL;
	uint32_t HV_Switch_State = 0;
    uint32_t HV_State[32] = {0};
	std::string HV_State_str = "0";
	std::string hv_filename = "hv_state_tmp.txt";
    hFile = fopen(hv_filename.c_str(),"rb");
    if(hFile == NULL){
        std::cout << "Could not open HV State file. Does it exist?" << std::endl;
        std::cout << "Error code: " << errno << std::endl;
    }else{
        fseek(hFile, 0, SEEK_END);
        int hSize = ftell(hFile);
        if(hSize != 128){
            std::cout << "HV State File Size Mismatch. Size: " << hSize << "\tSkipping this sample"<< std::endl;     
        }else{
            rewind(hFile);
            hBuffer = (uint32_t*) malloc (sizeof(uint32_t)*32);
            size_t hResult = fread(hBuffer, 1, hSize, hFile);
            for(int i=0; i<32; i++){
            	HV_State[i] = *(hBuffer + i);
				if(HV_State[i] == 2){
					HV_Switch_State = (HV_Switch_State | (1 << i));
				}else{
					std::cout << "Unknown HV State for SIAB#" << i << ": " << HV_State[i] << std::endl;
				}
            }
        }
        HV_State_str = Word_to_Str(HV_Switch_State);

        fclose(hFile);
        free(hBuffer);
    }
    system("sudo rm hv_state_tmp.txt &");
    usleep(1000);
    return HV_State_str;
}


std::string Get_Camera_Temp()
{
    std::cout << "Getting Camera Temperature" << std::endl;
    wqsiab.send(GET_ALL_SIAB_SIPM_TEMP);
    sleep(2);
    wqsiab.send(GET_ALL_SIAB_UC_TEMP);
    sleep(2);

    FILE *sFile = NULL;
    uint32_t *sBuffer = NULL;
    uint32_t SiPM_temp[32] = {0};
    std::string Temp_str = "0";
    std::string sipm_filename = "sipm_tmp.txt";
    sFile = fopen(sipm_filename.c_str(),"rb");
    if(sFile == NULL){
        std::cout << "Could not open SiPM Temp file. Does it exist?" << std::endl;
        std::cout << "Error code: " << errno << std::endl;
    }else{
        fseek(sFile, 0, SEEK_END);
        int sSize = ftell(sFile);
        if(sSize != 128){
            std::cout << "SiPM Temp File Size Mismatch. Size: " << sSize << "\tSkipping this sample"<< std::endl;       
        }else{
            rewind(sFile);
            sBuffer = (uint32_t*) malloc (sizeof(uint32_t)*32);
            size_t sResult = fread(sBuffer, 1, sSize, sFile);
            for(int i=0; i<32; i++){
                SiPM_temp[i] = *(sBuffer+i);
            }
        }
        fclose(sFile);
        free(sBuffer);
    }

    FILE *uFile = NULL;
    uint32_t *uBuffer = NULL;
    uint32_t UC_temp[32] = {0};
    std::string UC_temp_str = "0";
    std::string uc_filename = "uc_tmp.txt";
    uFile = fopen(uc_filename.c_str(),"rb");
    if(uFile == NULL){
        std::cout << "Could not open UC Temp file. Does it exist?" << std::endl;
        std::cout << "Error code: " << errno << std::endl;
    }else{
        fseek(uFile, 0, SEEK_END);
        int uSize = ftell(uFile);
        if(uSize != 128){
            std::cout << "UC Temp File Size Mismatch. Size: " << uSize << "\tSkipping this sample"<< std::endl;     
        }else{
            rewind(uFile);
            uBuffer = (uint32_t*) malloc (sizeof(uint32_t)*32);
            size_t uResult = fread(uBuffer, 1, uSize, uFile);
            for(int i=0; i<32; i++){
                UC_temp[i] = *(uBuffer+i);
            }
        }
        fclose(uFile);
        free(uBuffer);
    }

    for(int i=0; i<16; i++){
            std::string stemp_str = Byte_to_Str(SiPM_temp[SiPM_Temp_List[i]]/4);
            std::string utemp_str = Byte_to_Str(UC_temp[UC_Temp_List[i]]/4);
            Temp_str = Temp_str + stemp_str + utemp_str;
    }

    system("sudo rm sipm_tmp.txt &");
    usleep(1000);
    system("sudo rm uc_tmp.txt &");
    usleep(1000);
    return Temp_str;
}

void Get_State_MSG()
{
	if(EN_STATE_MSG)
	{
		std::string state_msg = "0";

		// This step adds the status of the processes to state message
		state_msg = "CT" + Byte_to_Str(SIAB_STATE) + Byte_to_Str(TRGB_STATE) + Byte_to_Str(COBO_STATE) + Byte_to_Str(LVPS_STATE) + Byte_to_Str(PDU_STATE) + "0";

		// This step adds the UNIX timestamp to the state message
		unsigned long int sec = time(NULL);
		std::string sec_str = Word_to_Str((uint32_t)sec);
		state_msg = state_msg + sec_str;

		// This step adds the EMON data to the state message
		int a,b = 0;
		Get_EMON(a,b);
		std::string em1 = Short_to_Str(a);
		std::string em2 = Short_to_Str(b);
		state_msg = state_msg + em1 + em2;

		state_msg = state_msg + Get_PDU_State();

		// This step adds the HK data to the state message
		uint32_t* hk_data;
		hk_data = Get_LVPS_data();
		uint32_t Music_PWR_State = 0;
		uint32_t SIAB_current[32] = {0};
		uint8_t HV_current[8] = {0};
		uint32_t HV_value[8] = {0};
		//uint32_t MEB_current[8] = {0};
		//uint32_t MEB_temp[8] = {0};
		std::string HV_Value_str;
		std::string HV_Current_str;
		for(int i=0; i<32; i++){
			SIAB_current[i] = *(hk_data+i);
			if(SIAB_current[i] > 4096000){
				Music_PWR_State = (Music_PWR_State | (1 << i));
			}
		}
		std::string Music_State = Word_to_Str(Music_PWR_State);

		for(int i=0; i<8; i++){
			HV_value[i] = *(hk_data+i+48);
			HV_Value_str = HV_Value_str + Word_to_Str(HV_value[i]);
			HV_current[i]  = (uint8_t)((*(hk_data+i+32))*2.441406E-06);
			HV_Current_str = HV_Current_str + Byte_to_Str(HV_current[i]);
			//MEB_current[i] = *(hk_data+i+40);	// will include later, if needed
			//MEB_temp[i]    = *(hk_data+i+56);	// will include later, if needed
		}

		state_msg += Music_State + Get_Camera_HV_State() + HV_Value_str + HV_Current_str + Get_Camera_Temp();
		wqcs.send(state_msg);
		std::cout << "STATE MESSAGE: " << state_msg << std::endl;
	}
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


////////////////////////////////////////////////////////
///////////////// Trigger Scan Commands ////////////////
////////////////////////////////////////////////////////

std::string Update_TB_GC_Scan(int Run_Num)
{
	std::string temp = fcutils::string::int_to_hex(Run_Num);
	std::string run_num_str = temp.substr(6, 4);
	std::string GC_CMD = "0A01" + run_num_str + "00000000";
	return GC_CMD;
}

std::string Update_TB_CC_Scan(int Run_Num)
{
	std::string temp = fcutils::string::int_to_hex(Run_Num);
	std::string run_num_str = temp.substr(6, 4);
	std::string CC_CMD = "0B01" + run_num_str + "00000000";
	return CC_CMD;
}

void Do_Trigger_Scan(std::string& msg)
{
	std::string NofSteps_str = msg.substr(0, 2);
	int NofSteps = stoi(NofSteps_str,0,16);

	std::string StepSize_str = msg.substr(2, 2);
	int StepSize = stoi(StepSize_str,0,16);

	std::string StepDuration_str = msg.substr(4, 2);
	int StepDuration = stoi(StepDuration_str,0,16);

	std::string SIAB_IDs_str = msg.substr(6, 8);

	std::string TScan_Get_GC = Update_TB_GC_Scan(TScanNo);	// This just updates the TB command with run number
	std::string TScan_Get_CC = Update_TB_CC_Scan(TScanNo);	// This just updates the TB command with run number
	wqtrgb.send(TB_Init);
	sleep(1);
	
	for(int i=0; i<NofSteps; i++)
	{
		int trg_threshold = 1 + i*StepSize;
		std::stringstream trg_stream;
		trg_stream << std::setfill('0') << std::setw(4) << std::hex << trg_threshold;
		std::string siab_cmd = SIAB_IDs_str + "2503" + trg_stream.str();
		wqsiab.send(siab_cmd);
		sleep(5);

		wqtrgb.send(TB_Prepare);
		sleep(1);
		wqtrgb.send(TB_Tscan_Config);
		sleep(2);
		wqtrgb.send(TB_Enable_Trigger);
		for (int j=0; j<(StepDuration/2); j++)
		{
			sleep(2);
			wqtrgb.send(TScan_Get_CC);
		}
		sleep(2);
		wqtrgb.send(TB_Disable_Trigger);
		sleep(1);
		wqtrgb.send(TScan_Get_GC);
		sleep(2);
		wqtrgb.send(TB_Clear_Counters);
		sleep(1);
	}
	std::cout << "Finished Trigger Scan Number:" << TScanNo << std::endl;
	TScanNo++;
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
///////////////// Control SW Commands //////////////////
////////////////////////////////////////////////////////

void Process_RC_Msg(const std::string& msg)
{
	std::string CMD_str = msg.substr(0, 2);
	std::string msg_payload = msg.substr(2, 14);
	int CMD_ID = stoi(CMD_str,0,16);
	//std::cout << "CMD_ID: " << CMD_ID << " was received." << std::endl;

	switch(CMD_ID)
	{
		case REBOOT_CPU_ID:
			Reset_CT_CPU();
			break;
		case SHUTDOWN_CPU_ID:
			Shutdown_CT_CPU();
			break;
		case Init_SIAB_ID:
			Init_SIAB();
			break;
		case Init_TRGB_ID:
			Init_TRGB();
			break;
		case Init_COBO_ID:
			Init_COBO();
			break;
		case Init_LVPS_ID:
			Init_LVPS();
			break;
		case Init_PDU_ID:
			Init_PDU();
			break;
		case Flush_FTP_ID:
			Flush_FTP();
			break;
		case Update_Runtime_ID:
			Update_Runtime(msg_payload);
			break;
		case RUN_START:
			Start_Run();
			break;
		case RUN_STOP:
			Stop_Run();
			break;
		case STATE_MSG_EN_ID:
			Enable_State_MSG();
			break;
		case STATE_MSG_DIS_ID:
			Disable_State_MSG();
			break;
		case GET_LOG_ID:
			Get_LOG(msg_payload);
			break;
		case DO_TRG_SCAN:
			Do_Trigger_Scan(msg_payload);
			break;
		case Kill_SIAB_ID:
			Kill_SIAB();
			break;
		case Kill_TRGB_ID:
			Kill_TRGB();
			break;
		case Kill_COBO_ID:
			Kill_COBO();
			break;
		case Kill_LVPS_ID:
			Kill_LVPS();
			break;
		case Kill_PDU_ID:
			Kill_PDU();
			break;
		case Kill_ALL_ID:
			Kill_All();
			break;
		default:
			break;
	}
}

void Process_Incoming_CMD(const std::string& msg)
{
	std::string qid_str = msg.substr(0, 2);
	std::string msg_payload = msg.substr(2, 16);
	int qid = stoi(qid_str,0,16);

	switch(qid)
	{
		case QID_RC:
			Process_RC_Msg(msg_payload);
			break;
		case QID_SIAB:
			wqsiab.send(msg_payload);
			break;
		case QID_TRGB:
			wqtrgb.send(msg_payload);
			break;
		case QID_COBO:
			wqcobo.send(msg_payload);
			break;
		case QID_LVPS:
			wqlvps.send(msg_payload);
			break;
		case QID_PDU:
			wqpdu.send(msg_payload);
			break;
		default:
			break;
	}
}

void read_mq_cs(sigval_t sig){
	int numRead{-1};
	// Retrive the posix message queue object from the sigval_t infos.
	posix_message_queue* mqp=static_cast<posix_message_queue*>(sig.sival_ptr);
	auto msg_size=mqp->get_msgsize();
	auto msg=new char[msg_size];
	auto mq_name=mqp->get_name();
	mqp->notify(read_mq_cs);
	while ((numRead = mqp->receive(msg, msg_size)) > 0){
		std::cout << " Queue: " << mq_name <<", message received: " << msg << std::endl;
	}
	Process_Incoming_CMD(msg);
	msg[0]='\0';
	delete [] msg;
}

void read_mq_rc(sigval_t sig){
	int numRead{-1};
	// Retrive the posix message queue object from the sigval_t infos.
	posix_message_queue* mqp=static_cast<posix_message_queue*>(sig.sival_ptr);
	auto msg_size=mqp->get_msgsize();
	auto msg=new char[msg_size];
	auto mq_name=mqp->get_name();
	mqp->notify(read_mq_rc);
	while ((numRead = mqp->receive(msg, msg_size)) > 0){
		std::cout << " Queue: " << mq_name <<", message received: " << msg << std::endl;
	}

	state_response = msg;
	delete [] msg;
}

void master_timer(std::function<void(void)> func, uint32_t interval)
{
	std::thread([func, interval]()
	{
		while(true)
		{
			auto x = std::chrono::steady_clock::now() + std::chrono::milliseconds(interval);
			func();
			std::this_thread::sleep_until(x);
		}
	}).detach();
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

int main()
{
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::cout << "Staring Master Control at: " << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << std::endl;

	rqcs.open();
	wqcs.open();
	rqrc.open();
	wqsiab.open();
	wqtrgb.open();
	wqcobo.open();
	wqlvps.open();
	wqpdu.open();

	rqcs.notify(read_mq_cs);
	rqrc.notify(read_mq_rc);

	master_timer(Check_Run_Status, run_interval*1000);
	master_timer(Get_State_MSG, state_msg_interval*1000);

	while(true)
	{
		usleep(100000);
	}
	
	return 0;
}
