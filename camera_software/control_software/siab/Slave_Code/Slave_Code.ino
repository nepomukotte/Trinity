#include <SPI.h>
#include <EEPROM.h>
#include "SoftSPI.h"
#include "Music_chip.h"
#include "AD7173.h"
#include <YetAnotherPcInt.h>
#include "device_specific.h"
#include "SMB_Commands.h"
#include "SMB_Slave.h"

#define SERIAL_BAUDRATE_HZ  115200
#define CONV_FACTOR         0.000149011612  // It converts ADC counts to mV
#define DEBUG_STATUS           false

int volatile tempofUC = 0;

uint32_t volatile adc_sample_sum[16];
uint16_t volatile adc_counter[16] = {0};
uint32_t volatile adc_avg[16] = {0};
uint32_t volatile adc_ref_offset[16] = {0};

bool volatile ADC_RUN = false;              // This flag will be set to true by users, whenever they want to start ADC interrupt routine.
bool volatile ADC_OFFSET_STATUS = false;    // This flag will be set to true, the offset of ADC is recorded. This flag will be checked in ADC_ready to turn off HV, if needed.

uint8_t  volatile music_pin_num;
uint8_t  volatile spi_addrs_bot[23];
uint8_t  volatile spi_addrs_top[23];
uint16_t volatile spi_data_bot[23];
uint16_t volatile spi_data_top[23];

const char *music_output[] = {"BCOUT0", "BCOUT1", "BCOUT2", "BCOUT3", "BCOUT4", "BCOUT5", "BCOUT6", "BCOUT7",
                              "TCOUT7", "TCOUT6", "TCOUT5", "TCOUT4", "TCOUT3", "TCOUT2", "TCOUT1", "TCOUT0"};

void setup() {

  delay(3000);

  Serial.begin(SERIAL_BAUDRATE_HZ);
  Serial.print("SIAB#");
  Serial.print(SMB_OWN_ADDRESS - 16);
  Serial.println(": Please Enter Your Command ...");
  
  pinMode(LED_BUILTIN, OUTPUT);

  pinMode(Music_3V3_Enable, OUTPUT);
  pinMode(Music_5V_Enable, OUTPUT);
  pinMode(musicChipCSBot, OUTPUT);
  pinMode(musicChipCSTop, OUTPUT);
  pinMode(music_top_rst, OUTPUT);
  pinMode(music_bot_rst, OUTPUT);
  pinMode(ADCCS, OUTPUT);
  pinMode(interruptPin, INPUT);
  pinMode(HV_Enable, OUTPUT);

  digitalWrite(Music_3V3_Enable, LOW);
  digitalWrite(Music_5V_Enable, LOW);
  digitalWrite(HV_Enable, LOW);

  digitalWrite(ADCCS, HIGH);
  digitalWrite(musicChipCSBot, HIGH);
  digitalWrite(musicChipCSTop, HIGH);
  digitalWrite(music_bot_rst, LOW);
  digitalWrite(music_top_rst, LOW);
  digitalWrite(music_bot_rst, HIGH);
  digitalWrite(music_top_rst, HIGH);
  digitalWrite(music_bot_rst, LOW);
  digitalWrite(music_top_rst, LOW);

  SMBusInit();
  SMBEnable();
  Setup_Timer();
}

void Setup_Timer()
{
  cli();                      // disable all interrupts
  TCCR1A = 0;
  TCCR1B = 0;
  TCNT1  = 0;
  OCR1A = 31250;              // compare match register 16MHz/256/1Hz
  TCCR1B |= (1 << WGM12);     // CTC mode
  TCCR1B |= (1 << CS12);      // 256 prescaler 
  TIMSK1 |= (1 << OCIE1A);    // enable timer compare interrupt
  sei();                      // enable all interrupts
}

ISR(TIMER1_COMPA_vect)
{
  if(ADC_RUN)
  {
    PcInt::detachInterrupt(interruptPin);
    for (int i=0; i<ENABLED_CH; i++)
    {
      adc_avg[i] = (adc_sample_sum[i] / adc_counter[i]);
      if(DEBUG_STATUS) {
        Serial.print("SiPM#");
        Serial.print(i+1);
        Serial.print(": ");
        Serial.print(adc_avg[i]*CONV_FACTOR, 2);
        Serial.print("mV [");
        Serial.print(music_output[i]);
        Serial.println("]");
        }
      adc_counter[i] = 0;
      adc_sample_sum[i] = 0;
    }
    if(DEBUG_STATUS) {Serial.println("");}
    PcInt::attachInterrupt(interruptPin, ADC_ready, (const char*) 0 , FALLING);
  }
}

void GetTemp(void)
{
  uint32_t wADC;
  ADMUX = (_BV(REFS1) | _BV(REFS0) | _BV(MUX3));  // Set the voltage reference and mux.
  ADCSRA |= _BV(ADEN);                            // enable the ADC
  delay(50);                                      // wait for voltages to become stable.
  ADCSRA |= _BV(ADSC);                            // Start the ADC
  while (bit_is_set(ADCSRA,ADSC));                // Detect end-of-conversion
  wADC = ADCW;                                    // Reading register "ADCW" takes care of how to read ADCL and ADCH.
  tempofUC = (int)wADC;
  //tempofUC = (wADC - 324.31 ) / 1.22;             // The offset of 324.31 could be wrong. It is just an indication.
}

void loop() {
  GetTemp();
  delay(1000);
}
