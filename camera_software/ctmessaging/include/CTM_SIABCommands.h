#ifndef SIAB_H
#define SIAB_H

#define NSIABCOMMANDS 38
#define NMUSICOPTS 3

#define CMD_RESET               		    0x01
#define CMD_TURN_ON_LED         		    0x02
#define CMD_TURN_OFF_LED        		    0x03
#define CMD_TURN_ON_3V3         		    0x04
#define CMD_TURN_OFF_3V3        		    0x05
#define CMD_TURN_ON_5V          		    0x06
#define CMD_TURN_OFF_5V         		    0x07
#define CMD_TURN_ON_3V3_5V      		    0x08
#define CMD_TURN_OFF_3V3_5V     		    0x09
		  
#define CMD_READ_SIPM_TEMP      		    0x0A
#define CMD_READ_UC_TEMP        		    0X0B
		  
#define CMD_READ_EPRM_BYTE      		    0X0C
#define CMD_READ_EPRM_WORD      		    0X0D
#define CMD_WRITE_EPRM_BYTE     		    0X0E
#define CMD_WRITE_EPRM_WORD     		    0X0F

#define CMD_INIT_MUSIC                      0x10
#define CMD_LOAD_MUSIC_CONFIG_EPRM          0x11
#define CMD_PRINT_LOADED_MUSIC_CONFIG       0x12
#define CMD_WRITE_LOADED_CONFIG_TO_MUSIC    0x13
#define CMD_READ_FROM_MUSIC_REGISTER        0x14
#define CMD_WRITE_TO_MUSIC_REGISTER         0x15

#define CMD_SET_ADC_RATE         			0x16
#define CMD_INITIALIZE_ADC       			0x17
#define CMD_ENABLE_ADC_RUN       			0x18
#define CMD_DISABLE_ADC_RUN      			0x19
#define CMD_READ_SINGLE_CURRENT  			0x1A
#define CMD_READ_ALL_CURRENT     			0x1B
#define CMD_RECORD_ADC_OFFSET    			0x1C
#define CMD_READ_ADC_OFFSET      			0x1D

#define CMD_TURN_ON_HV           			0x1E
#define CMD_TURN_OFF_HV          			0x1F
#define CMD_ENABLE_HV_CTR        			0x20
#define CMD_CHECK_HV_STATUS      			0x21
#define CMD_SET_MAX_CURRENT		 			0x22

#define CMD_UPDATE_MUSIC_CONFIG	 			0x23
#define CMD_SET_BIAS_VOLTAGE	 			0x24
#define CMD_SET_TRG_THRESHOLD				0x25

#define CMD_PING_SIAB						0x30

typedef struct{
	const char *names[38]; 
	int values[38]; 
} siabDict;

typedef struct{
	const char *names[3];
	int values[3];
} musicDict;

siabDict siabCommands = { .names = {"reset\0",
					"powerLEDOn\0",
					"powerLEDOff\0",
					"power3v3On\0",
					"power3v3Off\0",
					"power5vOn\0",
					"power5vOff\0",
					"powerOn\0",
					"powerOff\0",

					"sipmTemp\0",
					"ucTemp\0",

					"readEPRMByte\0",
					"readEPRMWord\0",
					"writeEPRMByte\0",
					"writeEPRMWord\0",

					"MusicInit\0",
					"MusicConfLoadEPRM\0",
					"MusicConfPrint\0",
					"MusicConfWrite\0",
					"MusicRegRead\0",
					"MusicRegWrite\0",
					"ADCsetRate\0",
					"ADCInitialize\0",
					"ADCRunEnable\0",
					"ADCRunDisable\0",
					"ADCReadCurrentSingle\0",
					"ADCReadCurrentAll\0",
					"ADCRecordOffset\0",
					"ADCReadOffset\0",
					"powerHVON\0",
					"powerHVOFF\0",
					"HVControlEnable\0",
					"HVControlStatus\0",
					"HVControlSetMaxCurrent\0",
					"MusicConfUpdate\0",
					"loadBiasVoltage\0",
					"MusicConfTrigTh\0",
					"ping\0"},
				.values = {CMD_RESET,               		
					CMD_TURN_ON_LED,         		
					CMD_TURN_OFF_LED,        		
					CMD_TURN_ON_3V3,         		
					CMD_TURN_OFF_3V3,        		
					CMD_TURN_ON_5V,          		
					CMD_TURN_OFF_5V,         		
					CMD_TURN_ON_3V3_5V,      		
					CMD_TURN_OFF_3V3_5V,     		
					  
					CMD_READ_SIPM_TEMP,      		
					CMD_READ_UC_TEMP,        		
					  
					CMD_READ_EPRM_BYTE,      		
					CMD_READ_EPRM_WORD,      		
					CMD_WRITE_EPRM_BYTE,     		
					CMD_WRITE_EPRM_WORD,

					CMD_INIT_MUSIC,                  
					CMD_LOAD_MUSIC_CONFIG_EPRM,      
					CMD_PRINT_LOADED_MUSIC_CONFIG,   
					CMD_WRITE_LOADED_CONFIG_TO_MUSIC,
					CMD_READ_FROM_MUSIC_REGISTER,    
					CMD_WRITE_TO_MUSIC_REGISTER,     
					CMD_SET_ADC_RATE,         		
					CMD_INITIALIZE_ADC,       		
					CMD_ENABLE_ADC_RUN,       		
					CMD_DISABLE_ADC_RUN,      		
					CMD_READ_SINGLE_CURRENT,  		
					CMD_READ_ALL_CURRENT,     		
					CMD_RECORD_ADC_OFFSET,    		
					CMD_READ_ADC_OFFSET,      		
					CMD_TURN_ON_HV,           		
					CMD_TURN_OFF_HV,          		
					CMD_ENABLE_HV_CTR,        		
					CMD_CHECK_HV_STATUS,      		
					CMD_SET_MAX_CURRENT,		 		
					CMD_UPDATE_MUSIC_CONFIG,	 		
					CMD_SET_BIAS_VOLTAGE,
					CMD_SET_TRG_THRESHOLD,
					CMD_PING_SIAB	 		
					}};

musicDict musicIDs = { .names = {"bottom\0",
					"top\0",
					"both\0"},
					.values = {1,2,3}};
#endif