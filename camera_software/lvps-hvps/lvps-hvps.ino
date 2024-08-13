#include <SPI.h>
#include "AD7173.h"
#include "AD5686R.h"
#include <YetAnotherPcInt.h>

// ***** ADC Initialization ***** //
#define NofADC 4
#define ENABLED_CH 16
#define interruptPin 12
#define max_adc_value 16777215
#define conv_factor 0.00000596046     //this is the result of (100.0 mA/max_adc_value).

const int log_data = 1;
volatile int ADCFlag = 0;
volatile int sample_counter = 0;
volatile double adc_vals[NofADC][16];
const char *HK_list[] = {"BP0", "BP1", "HK1", "HK2"}; 
// ****************************** //

// ***** DAC Initialization ***** //
#define CSDAC 10
float HV1_value = 44.0;
float HV2_value = 44.0;
float HV3_value = 44.0;
float HV4_value = 44.0;
unsigned short voltage[4] = {0, 0, 0, 0};
// ****************************** //

void setup()
{
  delay(100);
  Serial.begin(115200);
  pinMode(CSDAC, OUTPUT);
  digitalWrite(CSDAC, HIGH);
  delay(100);

  pinMode(interruptPin, INPUT);
  for (int i=0; i<NofADC; i++){
    pinMode(CSHK[i], OUTPUT);
    digitalWrite(CSHK[i], HIGH);
    delay(500);
  }

  AD5686R_Init();
  AD5686R_Reset();
  AD5686R_WriteFunction(AD5686R_CMD_SET_DCEN, AD5686R_ADDR_DAC_ALL, AD5686R_DCEN_ENABLE);
  AD5686R_InternalVoltageReference(AD5686R_INT_REF_ON);
  AD5686R_SetVoltage(AD5686R_ADDR_DAC_A, HV1_value);
  AD5686R_SetVoltage(AD5686R_ADDR_DAC_B, HV2_value);
  AD5686R_SetVoltage(AD5686R_ADDR_DAC_C, HV3_value);
  AD5686R_SetVoltage(AD5686R_ADDR_DAC_D, HV4_value);

  AD7173.init();
  for (int i=0; i<NofADC; i++){
    AD7173.init_ADC(i+1);
    delay(1000);
  }
  
  PcInt::attachInterrupt(interruptPin, ADC_ready_int, (const char*) 0 , FALLING);
  delay(1000);
}

void ADC_ready_int(const char* message, bool pinstate)
{
  cli();
  if (sample_counter < 16) {
  byte result[4];
  AD7173.get_data(result);
  if (result[3]&0x70) {return;}       //check for any error bits
  uint8_t channel = result[3]&0x0f;   //channel is last 4 bits
  adc_vals[ADCFlag][channel] = (((uint32_t)result[0] << 16) | ((uint32_t)result[1] << 8) | ((uint32_t)result[2]));
  sample_counter++;
  }
  sei();
}

void get_current() {
  for (int i=0; i<NofADC; i++){
    ADCFlag = i;
    sample_counter = 0;
    digitalWrite(CSHK[i], LOW);
    delay(32);
    digitalWrite(CSHK[i], HIGH);
    }
  if(log_data) {
    Serial.print("\t");
    for (int i=0; i<ENABLED_CH; i++){
      Serial.print(i+1);
      Serial.print("\t");
    }
    Serial.println("");
    for (int i=0; i<NofADC; i++){
      Serial.print(HK_list[i]);
      Serial.print("\t");
      for (int j=0; j<ENABLED_CH; j++){
        Serial.print((adc_vals[i][j]*conv_factor));
        Serial.print("\t");
        }
      Serial.println("");
    }
    Serial.println("");
    memset(adc_vals, 0, sizeof(adc_vals));
  }
}

void loop()
{
  get_current();
  delay(1000);
}
