/******************************************************************************/
/***************************** Include Files **********************************/
/******************************************************************************/
#include "AD5686R.h"		// AD5686R definitions.

/******************************************************************************/
/************************ Variables Definitions *******************************/
/******************************************************************************/
unsigned char currentPowerRegValue = 0; 

/******************************************************************************/
/************************ Functions Definitions *******************************/
/******************************************************************************/

/***************************************************************************//**
 * @brief Initializes the device over SPI interface.
 *
 * @return status - Result of the initialization procedure.
 *					Example: 0x0 - SPI peripheral was not initialized.
 *				  			   0x1 - SPI peripheral is initialized.
*******************************************************************************/
unsigned char AD5686R_Init()
{
  unsigned char status = -1;
  status = SPI_Init();
  if (status == 0)
  {
    Serial.println("SPI Initialization was successful.");
  }
  else
  {
    Serial.println("Could not Initialize the SPI.");
  }
  return status;
}

/***************************************************************************//**
 * @brief Puts the device in a specific power mode.
 *
 * @param channel - Channel option.
 *                  Example: AD5686R_ADDR_DAC_A
 *                           AD5686R_ADDR_DAC_B
 *                           AD5686R_ADDR_DAC_C
 *                           AD5686R_ADDR_DAC_D
 * @param pwrMode - Power mode of the device.
 *                  Example: AD5686R_PD_NORMAL   - Normal operation
 *                           AD5686R_PD_1K       - 1 kOhm to GND
 *                           AD5686R_PD_100K     - 100 kOhm to GND
 *                           AD5686R_PD_3STATE   - Three-state
 *                               
 *
 * @return none.
*******************************************************************************/
void AD5686R_PowerMode(unsigned char channel, unsigned char pwrMode)
{    
    switch(channel)
    {
        case AD5686R_ADDR_DAC_A:
            currentPowerRegValue &= ~AD5686R_PWR_PDA(0x3);   // Clear power bits
            Serial.println(AD5686R_PWR_PDA(0x3));
            Serial.println(~AD5686R_PWR_PDA(0x3));
            Serial.println(currentPowerRegValue);
            currentPowerRegValue |= AD5686R_PWR_PDA(pwrMode);
            Serial.println(currentPowerRegValue);
            break;
        case AD5686R_ADDR_DAC_B:
            currentPowerRegValue &= ~AD5686R_PWR_PDB(0x3);
            currentPowerRegValue |= AD5686R_PWR_PDB(pwrMode);
            break;
        case AD5686R_ADDR_DAC_C:
            currentPowerRegValue &= ~AD5686R_PWR_PDC(0x3);
            currentPowerRegValue |= AD5686R_PWR_PDC(pwrMode);
            break;
        case AD5686R_ADDR_DAC_D:
            currentPowerRegValue &= ~AD5686R_PWR_PDD(0x3);
            currentPowerRegValue |= AD5686R_PWR_PDD(pwrMode);
            break;
    }
    
    AD5686R_SetInputRegister(AD5686R_CMD(AD5686R_CMD_POWERMODE) | 
                             currentPowerRegValue);
}

/***************************************************************************//**
 * @brief Resets the device(clears the outputs to zero scale).
 * 
 * This function was supposed to control the output of the DAC after the reset
 * or at startup based on the status of RSTSEL pin (RSTSEL->HIGH sets the DAC
 * output at mid scale or 2.5V and RSTSEL->LOW sets the DAC output at zero volts),
 * but because in our design we had the RSTSEL pin tied to ground, output will be
 * always zero at startup or after software reset.
 * Here we only transmit the software reset command (0110) to device.
 *
 * @return none.
*******************************************************************************/
void AD5686R_Reset()
{
  AD5686R_SetInputRegister(AD5686R_CMD(AD5686R_CMD_SOFT_RESET));
}

/***************************************************************************//**
 * @brief Writes a 24-bit data-word to the Input Register of the device.
 *
 * @param registerValue - Value of the register.
 *
 * @return none.
*******************************************************************************/
void AD5686R_SetInputRegister(unsigned long registerValue)
{
    unsigned char registerWord[3] = {0, 0, 0};
    unsigned char* dataPointer    = (unsigned char*)&registerValue;
    
    registerWord[0] = dataPointer[2];
    registerWord[1] = dataPointer[1];
    registerWord[2] = dataPointer[0];
    
    SPI_Write(registerWord);
}

/***************************************************************************//**
 * @brief Select internal or external voltage reference.
 *
 * @param vRefMode - Voltage reference option.
 *                   Example: AD5686R_INT_REF_ON  - Internal reference is used.
 *                            AD5686R_INT_REF_OFF - External reference is used.
 *
 * @return none.
*******************************************************************************/
void AD5686R_InternalVoltageReference(unsigned char vRefMode)
{
  AD5686R_SetInputRegister(AD5686R_CMD(AD5686R_CMD_INT_REF_SETUP) | vRefMode);
}

/***************************************************************************//**
 * @brief Write data to the Input Register or to DAC Register of a channel.
 *
 * @param writeCommand - The write command option.
 *                       Example: 
*                        AD5686R_CMD_WR_INPUT_N    - Write to Input Register n.
 *                       AD5686R_CMD_WR_UPDT_DAC_N - Write to and update DAC n.
 * @param channel - Channel option.
 *                  Example: AD5686R_ADDR_DAC_A
 *                           AD5686R_ADDR_DAC_B
 *                           AD5686R_ADDR_DAC_C
 *                           AD5686R_ADDR_DAC_D
 * @param data -  Data value to write.
 *
 * @return none.
*******************************************************************************/
void AD5686R_WriteFunction(unsigned char writeCommand, unsigned char channel, unsigned short data)
{
    AD5686R_SetInputRegister(AD5686R_CMD(writeCommand) | AD5686R_ADDR(channel) | (AD5686R_DATA_BITS(data)));
}

/***************************************************************************//**
 * @brief Reads back the binary value written to one of the channels.
 *
 * @param dacChannelAddr - Channel address.
 *                         Example: AD5686R_ADDR_DAC_A
 *                                  AD5686R_ADDR_DAC_B
 *                                  AD5686R_ADDR_DAC_C
 *                                  AD5686R_ADDR_DAC_D 
 *
 * @return 16-bit value of the selected channel.
*******************************************************************************/
unsigned short AD5686R_ReadBack(unsigned char dacChannelAddr)
{
  unsigned long channelValue = 0;
  unsigned char rxBuffer[3]  = {0, 0, 0};
  
  AD5686R_SetInputRegister(AD5686R_CMD(AD5686R_CMD_SET_READBACK) | AD5686R_ADDR(dacChannelAddr));
  SPI_Read(rxBuffer);
  channelValue = (((long)rxBuffer[1] << 8) | rxBuffer[2]);    
  return channelValue;
}

/***************************************************************************//**
 * @brief Selects the output voltage of the selected channel.
 *
 * @param channel - Channel option.
 *                  Example: AD5686R_ADDR_DAC_A
 *                           AD5686R_ADDR_DAC_B
 *                           AD5686R_ADDR_DAC_C
 *                           AD5686R_ADDR_DAC_D
 * @param outputVoltage - Output voltage value.
 * 
 * Please note that voltage reference can be adjusted as well based on the status
 * of GAIN pin, however since we have tied the GAIN pin to HIGH, the GAIN is set 
 * to 2 which results in a 5.0V reference voltage.
 *
 * @return The actual value of the output voltage.
*******************************************************************************/
float AD5686R_SetVoltage(unsigned char channel, float outputVoltage)
{
  float maxVoltage = 65.573;;
  float actualVoltage = 0;
  unsigned short binaryValue = 0;
  
  binaryValue = (unsigned short)(outputVoltage * (65536) / maxVoltage);
  AD5686R_WriteFunction(AD5686R_CMD_WR_UPDT_DAC_N, channel, binaryValue);
  actualVoltage = (maxVoltage * binaryValue) / 65536;
  
  return actualVoltage;
}
