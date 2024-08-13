// ***************************************************** //
// ********** Trigger Board Memory Map ***************** //
// ***************************************************** //

// ***************************************************** //
// ********** Useful Constant Parameters *************** //
// ***************************************************** //
#define All_Bits_Low				0x00000000
#define All_Bits_High				0xFFFFFFFF

#define Trig_Rate_1Hz				0x05F5E100
#define Trig_Rate_10Hz				0x00989680
#define Trig_Rate_100Hz				0x000F4240
#define Trig_Rate_1KHz				0x000186A0
#define Trig_Rate_10KHz				0x00002710
#define Trig_Rate_100KHz			0x000003E8
#define Trig_Rate_1MHz				0x00000064
#define Trig_Rate_10MHz				0x0000000A

#define Ten_ns_Interval				0x00000001
#define Hundred_ns_Interval			0x0000000A
#define One_us_Interval				0x00000064
#define Ten_us_Interval				0x000003E8
#define Hundred_us_Interval			0x00002710
#define One_ms_Interval				0x000186A0
#define Ten_ms_Interval				0x000F4240
#define Hundred_ms_Interval			0x00989680
#define One_Second_Interval			0x05F5E100
#define Ten_Second_Interval			0x3B9ACA00
#define Thirty_Second_Interval		0xB2D05E00

#define Dis_All_Triggers			0x00000000
#define En_All_Triggers				0xFFFFFFFF
#define En_BiFocal_Trig_only		0x00000001
#define En_Disc_Test_Trig_only		0x00000002
#define En_Int_Trig_only			0x00000004
#define En_Ext_Trig_only			0x00000008
#define En_GPS_Trig_only			0x00000010
#define En_BiF_And_Int_only			0x00000005

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

// ########## To be updated .. ######################### //


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
#define Overflow_1							0x200C
#define Overflow_2							0x2010
#define Overflow_3							0x2014

// ########## Counter Value Registers ################## //
#define Counters_First_Addr					0x2020
#define Clock_Counter_0						0x2020
#define Clock_Counter_1						0x2024
#define Event_Counter						0x2028

// ########## Trigger Component Counters ############### //
#define Bifocal_Counter_Addr				0x202C
#define Discriminator_Counter_Addr			0x2030
#define Internal_Counter_Addr				0x2034
#define External_Counter_Addr				0x2038
#define GPS_Counter_Addr					0x203C
#define LED_Counter_Addr					0x2040

// ########## Trigger Component Live Counters ########## //
#define Bifocal_Live_Counter_Addr			0x2044
#define Discriminator_Live_Counter_Addr		0x2048
#define Internal_Live_Counter_Addr			0x204C
#define External_Live_Counter_Addr			0x2050
#define GPS_Live_Counter_Addr				0x2054
#define LED_Live_Counter_Addr				0x2058

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