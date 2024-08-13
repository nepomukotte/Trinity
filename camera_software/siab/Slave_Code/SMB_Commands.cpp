#include <Arduino.h>
#include "SMB_Commands.h"
#include "SMB_Slave.h"
#include "Music_chip.h"
#include "AD7173.h"
#include <EEPROM.h>

AD7173Class AD7173;
MiniMusicBoard music;

// This flag will be set true, whenever the system is thermally stablized,
// and user has recorded a valid number in adc_ref_offset. In the ADC_ready
// function, the samples will be checked with adc_ref_offset values, and HV
// will be turned off, if they are larger than adc_ref_offset + 2mV.
bool HV_control_flag = false;

extern uint16_t tempofUC;

extern bool ADC_RUN;
extern bool ADC_OFFSET_STATUS;

data_rate_t ADC_Rate = SPS_2597;
const double adc_to_uamp = 38.8;
const uint16_t max_current = 400;
volatile double max_adc_count = adc_to_uamp * max_current;
uint32_t adc_sample_single[16];
extern uint32_t adc_avg[16];
extern uint32_t adc_ref_offset[16];
extern uint16_t adc_counter[16];
extern uint32_t adc_sample_sum[16];

extern uint8_t  music_pin_num;
extern uint8_t  spi_addrs_bot[23];
extern uint8_t  spi_addrs_top[23];
extern uint16_t spi_data_bot[23];
extern uint16_t spi_data_top[23];

void ProcessReceiveByte(SMBData *smb)
{
    smb->txBuffer[0] = ~PIND;
    smb->txLength = 1;
}

void ProcessMessage(SMBData *smb)
{
    if (smb->state == SMB_STATE_WRITE_REQUESTED)
    {
        switch (smb->rxBuffer[0])
        {
            case CMD_RESET:
                Reset_Atmega328P(smb);
                break;
            case CMD_TURN_ON_LED:
                Turn_on_LED(smb);
                break;
            case CMD_TURN_OFF_LED:
                Turn_off_LED(smb);
                break;
            case CMD_TURN_ON_3V3:
                Turn_on_3V3(smb);
                break;
            case CMD_TURN_OFF_3V3:
                Turn_off_3V3(smb);
                break;
            case CMD_TURN_ON_5V:
                Turn_on_5V(smb);
                break;
            case CMD_TURN_OFF_5V:
                Turn_off_5V(smb);
                break;
            case CMD_TURN_ON_3V3_5V:
                Turn_on_3V3_5V(smb);
                break;
            case CMD_TURN_OFF_3V3_5V:
                Turn_off_3V3_5V(smb);
                break;
            case CMD_READ_SIPM_TEMP:
                Read_SiPM_Temp(smb);
                break;
            case CMD_READ_UC_TEMP:
                Read_UC_Temp(smb);
                break;
            case CMD_READ_EPRM_BYTE:
                Read_EPRM_Byte(smb);
                break;
            case CMD_READ_EPRM_WORD:
                Read_EPRM_Word(smb);
                break;
            case CMD_WRITE_EPRM_BYTE:
                Write_EPRM_Byte(smb);
                break;
            case CMD_WRITE_EPRM_WORD:
                Write_EPRM_Word(smb);
                break;
            case CMD_INIT_MUSIC:
                Init_Music(smb);
                break;
            case CMD_LOAD_MUSIC_CONFIG_EPRM:
                Load_Music_Config_EPRM(smb);
                break;
            case CMD_PRINT_LOADED_MUSIC_CONFIG:
                Print_Loaded_Music_Config(smb);
                break;
            case CMD_WRITE_LOADED_CONFIG_TO_MUSIC:
                Write_Loaded_Config_to_Music(smb);
                break;
            case CMD_READ_FROM_MUSIC_REGISTER:
                Read_from_Music_Register(smb);
                break;
            case CMD_WRITE_TO_MUSIC_REGISTER:
                Write_to_Music_Register(smb);
                break;
            case CMD_SET_ADC_RATE:
                Set_ADC_Rate(smb);
                break;
            case CMD_INITIALIZE_ADC:
                Initialize_ADC(smb);
                break;
            case CMD_ENABLE_ADC_RUN:
                Enable_ADC_Run(smb);
                break;
            case CMD_DISABLE_ADC_RUN:
                Disable_ADC_Run(smb);
                break;
            case CMD_READ_SINGLE_CURRENT:
                Read_Single_Current(smb);
                break;
            case CMD_READ_ALL_CURRENT:
                Read_All_Current(smb);
                break;
            case CMD_RECORD_ADC_OFFSET:
                Record_ADC_Offset(smb);
                break;
            case CMD_READ_ADC_OFFSET:
                Read_ADC_Offset(smb);
                break;
            case CMD_TURN_ON_HV:
                Turn_on_HV(smb);
                break;
            case CMD_TURN_OFF_HV:
                Turn_off_HV(smb);
                break;
            case CMD_ENABLE_HV_CTR:
                Enable_HV_CTR(smb);
                break;
            case CMD_CHECK_HV_STATUS:
                Check_HV_Status(smb);
                break;
            case CMD_SET_MAX_CURRENT:
                Set_Max_Current(smb);
                break;
            default:
                Undefined_CMD(smb);
                break;                
        }
    }
    else
    {
        smb->state = SMB_STATE_IDLE;
    }
}

#pragma inline=forced
void Reset_Atmega328P(SMBData *smb)
{
  __asm__ __volatile__ ("jmp 0x0000");
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Turn_on_LED(SMBData *smb)
{
  digitalWrite(LED_Pin, HIGH);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Turn_off_LED(SMBData *smb)
{
  digitalWrite(LED_Pin, LOW);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Turn_on_3V3(SMBData *smb)
{
  digitalWrite(Music_3V3_Enable, HIGH);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Turn_off_3V3(SMBData *smb)
{
  digitalWrite(Music_3V3_Enable, LOW);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Turn_on_5V(SMBData *smb)
{
  digitalWrite(Music_5V_Enable, HIGH);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Turn_off_5V(SMBData *smb)
{
  digitalWrite(Music_5V_Enable, LOW);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Turn_on_3V3_5V(SMBData *smb)
{
  digitalWrite(Music_3V3_Enable, HIGH);
  digitalWrite(Music_5V_Enable, HIGH);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Turn_off_3V3_5V(SMBData *smb)
{
  digitalWrite(Music_3V3_Enable, LOW);
  digitalWrite(Music_5V_Enable, LOW);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Read_SiPM_Temp(SMBData *smb)
{
  int temperature = 0;
  temperature = analogRead(Temp_Pin);
  smb->txBuffer[0] = (temperature & 0xFF);
  smb->txBuffer[1] = ((temperature >> 8) & 0xFF);
  smb->txLength = 2;
  smb->state = SMB_STATE_WRITE_READ_REQUESTED;
}

#pragma inline=forced
void Read_UC_Temp(SMBData *smb)
{
  smb->txBuffer[0] = (tempofUC & 0xFF);
  smb->txBuffer[1] = ((tempofUC >> 8) & 0xFF);
  smb->txLength = 2;
  smb->state = SMB_STATE_WRITE_READ_REQUESTED;
}

#pragma inline=forced
void Read_EPRM_Byte(SMBData *smb)
{
  uint8_t value = 0;
  uint16_t addr = 0;
  addr = (((smb->rxBuffer[2]) << 8) | (smb->rxBuffer[1]));
  value = EEPROM.read(addr);
  smb->txBuffer[0] = value;
  smb->txBuffer[1] = 0;       // this is set to zero, so no random number would be sent in high byte.
  smb->txLength = 2;
  smb->state = SMB_STATE_WRITE_READ_REQUESTED;
}

#pragma inline=forced
void Read_EPRM_Word(SMBData *smb)
{
  uint8_t value_low_byte = 0;
  uint8_t value_high_byte = 0;
  uint16_t addr = 0;
  addr = (((smb->rxBuffer[2]) << 8) | smb->rxBuffer[1]);
  value_high_byte = EEPROM.read(addr);
  value_low_byte = EEPROM.read(addr+1);
  smb->txBuffer[0] = value_low_byte;
  smb->txBuffer[1] = value_high_byte;
  smb->txLength = 2;
  smb->state = SMB_STATE_WRITE_READ_REQUESTED;
}

// For writing one byte, 2nd and 3rd buffer are address and 4th one is byte.
// This way we can write to whole EEPROM.
#pragma inline=forced
void Write_EPRM_Byte(SMBData *smb)
{
  uint8_t value = 0;
  uint16_t addr = 0;
  addr = (((smb->rxBuffer[3]) << 8) | (smb->rxBuffer[2]));
  value = smb->rxBuffer[4];
  EEPROM.write(addr, value);  
  smb->state = SMB_STATE_IDLE;
}

// For writing one word, 2nd buffer is address and 3rd and 4th are word.
// This way we can write only to first 256 address of EEPROM.
#pragma inline=forced
void Write_EPRM_Word(SMBData *smb)
{
  uint8_t value_low_byte = 0;
  uint8_t value_high_byte = 0;
  uint8_t addr = 0;
  addr = smb->rxBuffer[2];
  value_low_byte = smb->rxBuffer[3];
  value_high_byte = smb->rxBuffer[4];
  EEPROM.write(addr, value_high_byte);
  EEPROM.write(addr+1, value_low_byte);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Init_Music(SMBData *smb)
{
  music.init(DEBUG_MODE , SOFTSPI_CLOCK_DIVIDER);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Load_Music_Config_EPRM(SMBData *smb)
{
  music_pin_num = smb->rxBuffer[1];
  if (music_pin_num == 1) { music.read_music_config_eeprom (music_pin_num, spi_addrs_bot, spi_data_bot);}
  if (music_pin_num == 2) { music.read_music_config_eeprom (music_pin_num, spi_addrs_top, spi_data_top);}
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Print_Loaded_Music_Config(SMBData *smb)
{
  music_pin_num = smb->rxBuffer[1];
  if (music_pin_num == 1)
  {
    Serial.println("Bottom music:");
    Serial.print("Address\t");
    Serial.println("Data");
    for (int i=0; i < 23; i++)
    {
      Serial.print(spi_addrs_bot[i]);
      Serial.print("\t");
      Serial.println(spi_data_bot[i], HEX);
    }
    Serial.print("\n");
  }
  else if (music_pin_num == 2)
  {
    Serial.println("Top music:");
    Serial.print("Address\t");
    Serial.println("Data");
    for (int i=0; i < 23; i++)
    {
      Serial.print(spi_addrs_top[i]);
      Serial.print("\t");
      Serial.println(spi_data_top[i], HEX);
    }
    Serial.print("\n");
  }
  else
  {
    Serial.println("Wrong music pin is received.");
  }
  
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Write_Loaded_Config_to_Music(SMBData *smb)
{
  music_pin_num = smb->rxBuffer[1];
  if (music_pin_num == 1) { music.configure_music(musicChipCSBot, spi_addrs_bot, spi_data_bot);}
  if (music_pin_num == 2) { music.configure_music(musicChipCSTop, spi_addrs_top, spi_data_top);}
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Read_from_Music_Register(SMBData *smb)
{
  uint8_t music_reg_addr = 0;
  uint16_t music_reg_value = 0;
  music_pin_num = smb->rxBuffer[2];     // note that this is the high byte in the data of process call.
  music_reg_addr = smb->rxBuffer[1];    // note that this is the low byte in the data of process call.
  if (music_pin_num == 1) { music_reg_value = music.spi_read( musicChipCSBot , music_reg_addr , 0, 0);}
  if (music_pin_num == 2) { music_reg_value = music.spi_read( musicChipCSTop , music_reg_addr , 0, 0);}
  smb->txBuffer[0] = (music_reg_value & 0xFF);
  smb->txBuffer[1] = ((music_reg_value >> 8) & 0xFF);
  smb->txLength = 2;
  smb->state = SMB_STATE_WRITE_READ_REQUESTED;
}

#pragma inline=forced
void Write_to_Music_Register(SMBData *smb)
{
  uint8_t music_reg_addr = 0;
  uint8_t value_low_byte = 0;
  uint8_t value_high_byte = 0;
  uint16_t music_reg_value = 0;
  music_pin_num = (((smb->rxBuffer[2]) >> 6) & 0xFF);    // this is for extracting 2 bits on left 
  music_reg_addr = ((smb->rxBuffer[2]) & 0x3F);         // this is for extracting 6 bits on right
  value_low_byte = smb->rxBuffer[3];
  value_high_byte = smb->rxBuffer[4];
  music_reg_value = ((uint16_t)value_high_byte << 8 | (uint16_t)value_low_byte);
  if (music_pin_num == 1) { music.spi_reg_write( musicChipCSBot , music_reg_addr , music_reg_value, false);}
  if (music_pin_num == 2) { music.spi_reg_write( musicChipCSTop , music_reg_addr , music_reg_value, false);}
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Set_ADC_Rate(SMBData *smb)
{
  uint8_t ADC_Rate_Code = smb->rxBuffer[1];
  switch(ADC_Rate_Code)
  {
    case 1:
      ADC_Rate = SPS_1;
      break;
    case 2:
      ADC_Rate = SPS_16;
      break;
    case 3:
      ADC_Rate = SPS_100;
      break;
    case 4:
      ADC_Rate = SPS_381;
      break;
    case 5:
      ADC_Rate = SPS_1007;
      break;
    case 6:
      ADC_Rate = SPS_2597;
      break;    
  }

  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Initialize_ADC(SMBData *smb)
{
  AD7173.init();
  AD7173.sync();
  AD7173.reset();
  if (AD7173.is_valid_id())
  {
    AD7173.set_channel_config(CH0,  true, SETUP0, AIN15, REF_NEG);  // BCOUT0
    AD7173.set_channel_config(CH1,  true, SETUP0, AIN14, REF_NEG);  // BCOUT1
    AD7173.set_channel_config(CH2,  true, SETUP0, AIN13, REF_NEG);  // BCOUT2
    AD7173.set_channel_config(CH3,  true, SETUP0, AIN12, REF_NEG);  // BCOUT3
    AD7173.set_channel_config(CH4,  true, SETUP0, AIN11, REF_NEG);  // BCOUT4
    AD7173.set_channel_config(CH5,  true, SETUP0, AIN10, REF_NEG);  // BCOUT5
    AD7173.set_channel_config(CH6,  true, SETUP0, AIN9,  REF_NEG);  // BCOUT6
    AD7173.set_channel_config(CH7,  true, SETUP0, AIN8,  REF_NEG);  // BCOUT7
    AD7173.set_channel_config(CH8,  true, SETUP0, AIN3,  REF_NEG);  // TCOUT7
    AD7173.set_channel_config(CH9,  true, SETUP0, AIN2,  REF_NEG);  // TCOUT6
    AD7173.set_channel_config(CH10, true, SETUP0, AIN1,  REF_NEG);  // TCOUT5
    AD7173.set_channel_config(CH11, true, SETUP0, AIN16, REF_NEG);  // TCOUT4
    AD7173.set_channel_config(CH12, true, SETUP0, AIN7,  REF_NEG);  // TCOUT3
    AD7173.set_channel_config(CH13, true, SETUP0, AIN6,  REF_NEG);  // TCOUT2
    AD7173.set_channel_config(CH14, true, SETUP0, AIN5,  REF_NEG);  // TCOUT1
    AD7173.set_channel_config(CH15, true, SETUP0, AIN4,  REF_NEG);  // TCOUT0

    AD7173.set_setup_config(SETUP0, UNIPOLAR, AIN_BUF_DISABLE, REF_INT); // REF_PWR, REF_INT, REF_AIN
    AD7173.set_offset_config(OFFSET0, 0);
    AD7173.set_filter_config(FILTER0, ADC_Rate);
    AD7173.set_adc_mode_config(CONTINUOUS_CONVERSION_MODE, INTERNAL_CLOCK, REF_ENABLE);
    AD7173.set_interface_mode_config(false, true);
  }
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Enable_ADC_Run(SMBData *smb)
{
  ADC_RUN = true;
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Disable_ADC_Run(SMBData *smb)
{
  ADC_RUN = false;
  PcInt::detachInterrupt(interruptPin);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Read_Single_Current(SMBData *smb)
{
  uint8_t channel_number = 0;
  channel_number = smb->rxBuffer[1];    // note that for channel number, which ranges from 0 to 15, one byte is used which is the low byte of process call.
  smb->txBuffer[0] = ((adc_avg[channel_number] >> 8) & 0xFF);          // out of 24 bit of ADC, this is the 8 bit on the middle. The LSB is useless due to noise.
  smb->txBuffer[1] = ((adc_avg[channel_number] >> 16) & 0xFF);         // out of 24 bit of ADC, this is the 8 bit on the middle.
  smb->txLength = 2;
  smb->state = SMB_STATE_WRITE_READ_REQUESTED;
}

#pragma inline=forced
void Read_All_Current(SMBData *smb)
{
  // Doing nothing for now. Will add this feature later,
  // when I know how block read works in AVR316.
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Record_ADC_Offset(SMBData *smb)
{
  for(int i=0; i<ENABLED_CH; i++)
  {
    adc_ref_offset[i] = adc_avg[i];
  }
  ADC_OFFSET_STATUS = true;
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Read_ADC_Offset(SMBData *smb)
{
  uint8_t channel_number = 0;
  channel_number = smb->rxBuffer[1];    // note that for channel number, which ranges from 0 to 15, one byte is used which is the low byte of process call.
  smb->txBuffer[0] = ((adc_ref_offset[channel_number] >> 8) & 0xFF);          // out of 24 bit of ADC, this is the 8 bit on the middle. The LSB is useless due to noise.
  smb->txBuffer[1] = ((adc_ref_offset[channel_number] >> 16) & 0xFF);         // out of 24 bit of ADC, this is the 8 bit on the middle.
  smb->txLength = 2;
  smb->state = SMB_STATE_WRITE_READ_REQUESTED;
}

#pragma inline=forced
void Turn_on_HV(SMBData *smb)
{
  digitalWrite(HV_Enable, HIGH);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Turn_off_HV(SMBData *smb)
{
  digitalWrite(HV_Enable, LOW);
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Enable_HV_CTR(SMBData *smb)
{
  HV_control_flag = true;
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Check_HV_Status(SMBData *smb)
{
  int HV_Status = digitalRead(HV_Enable);
  if(HV_Status == LOW) { smb->txBuffer[0] = 0x01;}
  if(HV_Status == HIGH) { smb->txBuffer[0] = 0x02;}
  smb->txBuffer[1] = 0;
  smb->txLength = 2;
  smb->state = SMB_STATE_WRITE_READ_REQUESTED;
}

#pragma inline=forced
void Set_Max_Current(SMBData *smb)
{
  uint8_t threshold_lsb = 0;
  uint8_t threshold_msb = 0;
  threshold_msb = smb->rxBuffer[1];
  threshold_lsb = smb->rxBuffer[2];
  max_adc_count = adc_to_uamp * (((uint16_t)threshold_msb << 8) | ((uint16_t)threshold_lsb));
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void Undefined_CMD(SMBData *smb)
{
  smb->error = TRUE;
  smb->state = SMB_STATE_IDLE;
}

#pragma inline=forced
void ADC_ready(const char* message, bool pinstate)
{  
  cli();
  byte result[4];
  AD7173.get_data(result);
  if(result[3] & 0x70) {return;}            // check for any error bits
  uint8_t channel = (result[3] & 0x0F);     // channel is last 4 bits
  adc_sample_single[channel] = (((uint32_t)result[0] << 16) | ((uint32_t)result[1] << 8) | ((uint32_t)result[2]));
  if(HV_control_flag == true)
  {
    // This checks to make sure a valid threshold exists before shuting down HV.
    if((ADC_OFFSET_STATUS == true) & ( adc_sample_single[channel] > (adc_ref_offset[channel] + max_adc_count)))
    { digitalWrite(HV_Enable, LOW);}
  }
  adc_sample_sum[channel] += adc_sample_single[channel];
  adc_counter[channel]++;
  sei();
}
