#ifndef __AD5686R_H__
#define __AD5686R_H__

/******************************************************************************/
/***************************** Include Files **********************************/
/******************************************************************************/
#include "Communication.h"	// Communication definitions.

/******************************************************************************/
/******************************** AD5686R *************************************/
/******************************************************************************/

/* AD5686R Input Register */
#define AD5686R_CMD(x)                  ((unsigned long)(0x000F & (x)) << 20)
#define AD5686R_ADDR(x)                 ((unsigned long)(0x000F & (x)) << 16)
#define AD5686R_DATA_BITS(x)            ((unsigned long)(0xFFFF & (x)) <<  0)

/* AD5686R_CMD(x) options */
#define AD5686R_CMD_NOP                 0 // No operation.
#define AD5686R_CMD_WR_INPUT_N          1 // Write to Input Register n.
#define AD5686R_CMD_UPDATE_DAC_N        2 // Update DAC Register n.
#define AD5686R_CMD_WR_UPDT_DAC_N       3 // Write to and update DAC Channel n.
#define AD5686R_CMD_POWERMODE           4 // Power down/power up DAC.
#define AD5686R_CMD_LDAC_MASK           5 // Hardware LDAC mask register.
#define AD5686R_CMD_SOFT_RESET          6 // Software reset(power-on reset).
#define AD5686R_CMD_INT_REF_SETUP       7 // Internal reference setup register.
#define AD5686R_CMD_SET_DCEN            8 // Set up DCEN register.
#define AD5686R_CMD_SET_READBACK        9 // Set up readback register.

/* AD5686R_ADDR(x) options */
#define AD5686R_ADDR_DAC_A          1  // DAC A
#define AD5686R_ADDR_DAC_B          2  // DAC B
#define AD5686R_ADDR_DAC_C          4  // DAC C
#define AD5686R_ADDR_DAC_D          8  // DAC D
#define AD5686R_ADDR_DAC_A_B        3  // DAC A and DAC B
#define AD5686R_ADDR_DAC_ALL        15 // All DACs

/* Daisy-Chain Enable Register(DCEN) definition */
#define AD5686R_DCEN_DISABLE        0
#define AD5686R_DCEN_ENABLE         1

/* Power Setup Register definition */
#define AD5686R_PWR_PDD(x)      (((x) & 0x3) << 6)
#define AD5686R_PWR_PDC(x)      (((x) & 0x3) << 4)
#define AD5686R_PWR_PDB(x)      (((x) & 0x3) << 2)
#define AD5686R_PWR_PDA(x)      (((x) & 0x3) << 0)

/* AD5686R_PWR_PDn(x) options(n = A, B, C, D) */
#define AD5686R_PD_NORMAL       0 // Normal operation
#define AD5686R_PD_1K           1 // 1 kOhm to GND
#define AD5686R_PD_100K         2 // 100 kOhm to GND
#define AD5686R_PD_3STATE       3 // Three-state

/* LDAC Mask Register definition */
#define AD5686R_MSK_CH_A        1
#define AD5686R_MSK_CH_B        2
#define AD5686R_MSK_CH_C        4
#define AD5686R_MSK_CH_D        8

/* Reference Setup Register definition */
#define AD5686R_INT_REF_ON      0
#define AD5686R_INT_REF_OFF     1

/******************************************************************************/
/************************ Functions Declarations ******************************/
/******************************************************************************/

/*! Initializes the device. */
unsigned char AD5686R_Init();

/*! Resets the device(clears the outputs to zero scale). */
void AD5686R_Reset();

/*! Puts the device in a specific power mode. */
void AD5686R_PowerMode(unsigned char channel, unsigned char pwrMode);

/*! Select internal or external voltage reference. */
void AD5686R_InternalVoltageReference(unsigned char vRefMode);

/*!  Writes a 24-bit data-word to the Input Register of the device. */
void AD5686R_SetInputRegister(unsigned long registerValue);

/*! Write data to the Input Register or to DAC Register of a channel. */
void AD5686R_WriteFunction(unsigned char writeCommand, unsigned char channel, unsigned short data);

/*! Reads back the binary value written to one of the channels. */
unsigned short AD5686R_ReadBack(unsigned char dacChannelAddr);

/*! Selects the output voltage of the selected channel. */
float AD5686R_SetVoltage(unsigned char channel, float outputVoltage);

#endif // __AD5686R_H__
