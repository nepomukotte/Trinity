// ***************************************************** //
// ********** Trigger Board Memory Map ***************** //
// ***************************************************** //

// ***************************************************** //
// ********** Useful Constant Parameters *************** //
// ***************************************************** //
#define All_Bits_Low				0x00000000
#define All_Bits_High				0xFFFFFFFF
#define First_Bit_True				0x00000001

#define Trig_Rate_0Hz				0x00000000
#define Trig_Rate_1Hz				0x05F5E100
#define Trig_Rate_3Hz				0x01FCA055
#define Trig_Rate_5Hz				0x01312D00
#define Trig_Rate_10Hz				0x00989680
#define Trig_Rate_20Hz				0x004C4B40
#define Trig_Rate_30Hz				0x0032DCD5
#define Trig_Rate_50Hz				0x001E8480
#define Trig_Rate_100Hz				0x000F4240
#define Trig_Rate_200Hz				0x0007A120
#define Trig_Rate_500Hz				0x00030D40
#define Trig_Rate_1kHz				0x000186A0
#define Trig_Rate_2kHz				0x0000C350
#define Trig_Rate_10kHz				0x00002710
#define Trig_Rate_100kHz			0x000003E8
#define Trig_Rate_1MHz				0x00000064

#define En_BiFocal_Trig_only		0x00000001
#define En_Disc_Test_Trig_only		0x00000002
#define En_Int_Trig_only			0x00000004
#define En_Ext_Trig_only			0x00000008
#define En_GPS_Trig_only			0x00000010
#define En_All_Triggers				0x0000001F
#define Dis_All_Triggers			0x00000000

#define strech_length_10ns			0x00000001
#define strech_length_20ns			0x00000002
#define strech_length_50ns			0x00000005
#define strech_length_100ns			0x0000000A
#define strech_length_200ns			0x00000014
#define strech_length_300ns			0x0000001E
#define strech_length_500ns			0x00000032
#define strech_length_1us			0x00000064

#define Ten_ns_Interval				0x00000001
#define Hundred_ns_Interval			0x0000000A
#define One_us_Interval				0x00000064
#define Ten_us_Interval				0x000003E8
#define Hundred_us_Interval			0x00002710
#define One_ms_Interval				0x000186A0
#define Ten_ms_Interval				0x000F4240
#define Hundred_ms_Interval			0x00989680
#define One_Second_Interval			0x05F5E100
#define Two_Second_Interval			0x0BEBC200
#define Five_Second_Interval		0x1DCD6500
#define Ten_Second_Interval			0x3B9ACA00
#define Twenty_Second_Interval		0x77359400
#define Thirty_Second_Interval		0xB2D05E00

#define Deadtime_10ns				0x00000001
#define Deadtime_20ns				0x00000002
#define Deadtime_30ns				0x00000003
#define Deadtime_40ns				0x00000004
#define Deadtime_50ns				0x00000005
#define Deadtime_100ns				0x0000000A
#define Deadtime_150ns				0x0000000F
#define Deadtime_200ns				0x00000014
#define Deadtime_300ns				0x0000001E
#define Deadtime_400ns				0x00000028
#define Deadtime_500ns				0x00000032
#define Deadtime_600ns				0x0000003C
#define Deadtime_700ns				0x00000046
#define Deadtime_800ns				0x00000050
#define Deadtime_900ns				0x0000005A
#define Deadtime_1us				0x00000064

// ***************************************************** //
// ********** Start address of all registers *********** //
// ***************************************************** //
#define System_Reg_Start_Addr				0x0000
#define Trigger_Reg_Start_Addr				0x1000
#define Counter_Reg_Start_Addr				0x2000
#define Block_Memory_Control_Reg_Start_Addr	0x3000
#define Block_Memory_Access_Reg_Start_Addr	0x8000

// ***************************************************** //
// ********** System Registers ************************* //
// ***************************************************** //
#define Cookie					0x0100
#define Mesa_SW_Name_Low		0x0104
#define Mesa_SW_Name_High		0x0108
#define IDROM_Address_Offset	0x010C
#define CT_SW_Name_Low			0x0110
#define CT_SW_Name_High			0x0114
#define CT_SW_Date				0x0118
#define CT_SW_Version			0x011C
#define Test_Reg_Address		0x0F00

// ***************************************************** //
// ********** Trigger Registers ************************ //
// ***************************************************** //

// ########## Trigger Control and Status Registers ##### //
#define Disc_Stretch_Enable_0				0x1000
#define Disc_Stretch_Enable_1				0x1004
#define Disc_Stretch_Length					0x1008
#define Disc_Stretch_Restart				0x100C
#define Disc_Test_Trigger_Prescale			0x1010
#define Internal_Trigger_Prescale			0x1014
#define Internal_Trigger_Mode				0x1018
#define External_Trigger_Mode				0x101C
#define GPS_Trigger_Mode					0x1020
#define LED_Delay							0x1024
#define Transit_Busy_Length					0x1028
#define Enable_Busy							0x102C
#define Clear_Busy							0x1030
#define Save_Busy							0x1034
#define Read_Busy_Status					0x1038
#define Enable_Trigger_Types				0x103C
#define Enable_Triggering					0x1040
#define Disc_Deadtime						0x1044

// ########## Local Event Data Buffer ################## //
#define Event_Number						0x1080
#define Time_0								0x1084
#define Time_1								0x1088
#define Trigger_List						0x108C
#define Disc_Strech_1						0x1090
#define Disc_Strech_2						0x1094


// ***************************************************** //
// ********** Counter Registers ************************ //
// ***************************************************** //

// ########## Counter Control and Status Registers ##### //
#define Clear_Counters						0x2000
#define Save_Counters						0x2004
#define Rate_Counter_Period					0x2008
#define Overflow_0							0x200C
#define Overflow_1							0x2010
#define Overflow_2							0x2014

// ########## Counter Value Registers ################## //
#define Counters_First_Addr					0x2020
#define Clock_Counter_0						0x2020
#define Clock_Counter_1						0x2024
#define Event_Counter						0x2028

// ########## Trigger Comppnent Counters ############### //
#define Bi_Focal_Counter					0x202C
#define Discriminator_Counter				0x2030
#define Internal_Counter					0x2034
#define External_Counter					0x2038
#define GPS_Counter							0x203C
#define LED_Counter							0x2040

// ########## Trigger Comppnent Live Counters ########## //
#define Bi_Focal_Live_Counter				0x2044
#define Discriminator_Live_Counter			0x2048
#define Internal_Live_Counter				0x204C
#define External_Live_Counter				0x2050
#define GPS_Live_Counter					0x2054
#define LED_Live_Counter					0x2058

// ########## Busy Logic Counters Registers ############ //
#define TB_Busy_Counter_0					0x205C
#define TB_Busy_Counter_1					0x2060
#define CoBo_Busy_Counter_0					0x2064
#define CoBo_Busy_Counter_1					0x2068
#define Transit_Busy_Counter_0				0x206c
#define Transit_Busy_Counter_1				0x2070
#define Buffer_Busy_Counter_0				0x2074
#define Buffer_Busy_Counter_1				0x2078
#define Memory_Full_Counter_0				0x207c
#define Memory_Full_Counter_1				0x2080
#define Memory_Reading_Counter_0			0x2084
#define Memory_Reading_Counter_1			0x2088
#define Memory_Wrtiting_Counter_0			0x208C
#define Memory_Writing_Counter_1			0x2090

// ########## Discriminator Rate Counter Registers ##### //
#define IO_0_rate_Addr						0x2094
#define Find_Disc_Rate_Counter_Addr(IO_Number, address) do { \
		address = IO_0_rate_Addr + IO_Number * (4); \
    } while (0);

// ***************************************************** //
// ********** Data Storage Registers ******************* //
// ***************************************************** //

// ########## Block Memory Control and Status Registers  //
#define Block_Memory_Enable					0x3000
#define Memory_Start_Address				0x3004
#define Events_Written						0x3008
#define Readout_Start						0x300C
#define Memory_Block_Select					0x3010
#define Readout_Done						0x3014