#ifndef LVPS_H
#define LVPS_H

#define NLVPSCOMMANDS 12

#define RESET_uC		0x01
#define GET_uC_TEMP		0x02
#define INITILIZE_ADC		0x03
#define ENABLE_ADC_LOGGING	0x04
#define DISABLE_ADC_LOGGING	0x05
#define GET_ADC_VALUES		0X06
#define SET_LOOP_DELAY		0X07
#define INITILIZE_DAC		0X08
#define SET_FIXED_HV		0X09
#define SET_UPDATED_HV_VALUE	0X0A
#define SWITCH_BINARY_LOG	0X0B
#define SWITCH_ASCII_LOG	0X0C

typedef struct{
	const char *names[NLVPSCOMMANDS]; 
	int values[NLVPSCOMMANDS]; 
} lvpsDict;

lvpsDict lvpsCommands = { .names = {"resetMC\0",
					"readMCTemp\0",
					"initADC\0",
					"ADCLogOn\0",
					"ADCLogOff\0",
					"readADC\0",
					"setLoopDelay\0",
					"initDAC\0",
					"setFixedHV\0",
					"setUpdatedHV\0",
					"setLogBin\0",
					"setLogASCII\0"},
				.values = {RESET_uC,
					GET_uC_TEMP,
					INITILIZE_ADC,
					ENABLE_ADC_LOGGING,
					DISABLE_ADC_LOGGING,
					GET_ADC_VALUES,
					SET_LOOP_DELAY,
					INITILIZE_DAC,
					SET_FIXED_HV,
					SET_UPDATED_HV_VALUE,
					SWITCH_BINARY_LOG,
					SWITCH_ASCII_LOG	
					}};
#endif