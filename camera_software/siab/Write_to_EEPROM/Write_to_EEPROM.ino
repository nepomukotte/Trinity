#include <EEPROM.h>

uint8_t data_low;
uint8_t data_high;
uint8_t eeprm_start_addr_bot = 0;
uint8_t eeprm_start_addr_top = 69;

uint8_t spi_addrs[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,32,33,34,35,36,37,38};

uint16_t spi_data_bot[] = {0x4FF , 0x4FF , 0x4FF , 0x4FF , 0x4FF , 0x4FF , 0x4FF , 0x4FF ,
                           0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7,
                           0x9090, 0x5081, 0x799E, 0xF3FF, 0x9BEE, 0x03DC, 0x08CB };

uint16_t spi_data_top[] = {0x4FF , 0x4FF , 0x4FF , 0x4FF , 0x4FF , 0x4FF , 0x4FF , 0x4FF ,
                           0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7, 0xF3C7,
                           0x9090, 0x5081, 0x799E, 0xF3FF, 0x9BEE, 0x03DC, 0x08CB };

void setup() {
  
  Serial.begin(115200);
  
  for (int i = 0 ; i < EEPROM.length() ; i++) {
    EEPROM.write(i, 0);
  }
  Serial.println("EEPROM is cleared.");

  for (int i=0; i < 23; i++)
  {
    data_low = (spi_data_bot[i] & 0xFF);
    data_high = ((spi_data_bot[i] >> 8) & 0xFF);
    EEPROM.write(eeprm_start_addr_bot, spi_addrs[i]);
    EEPROM.write(eeprm_start_addr_bot + 1, data_high);
    EEPROM.write(eeprm_start_addr_bot + 2, data_low);
    eeprm_start_addr_bot += 3;
  }
  Serial.println("Bottom Music data is written to EEPROM.");


  for (int i=0; i < 23; i++)
  {
    data_low = (spi_data_top[i] & 0xFF);
    data_high = ((spi_data_top[i] >> 8) & 0xFF);
    EEPROM.write(eeprm_start_addr_top, spi_addrs[i]);
    EEPROM.write(eeprm_start_addr_top + 1, data_high);
    EEPROM.write(eeprm_start_addr_top + 2, data_low);
    eeprm_start_addr_top += 3;
  }
  Serial.println("Top Music data is written to EEPROM.");

  for (int i=0; i < 138; i++)
  {
    Serial.println(EEPROM.read(i), HEX);
  }

}

void loop() {

}
