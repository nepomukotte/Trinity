#include <EEPROM.h>
#include "Low_Power.h"
#include "SoftSPI.h"

#include "Music_chip.h"

//Global pins: << This needs to be removed. From now on , we use A0 and A1 to reset music chips.
#define ASIC_RST     8

//Expansion Pin numbers:
#define EXP_2 PORTD2
#define EXP_3 PORTD3
#define EXP_4 PORTD4
#define EXP_5 PORTD5
#define EXP_6 PORTD6
#define EXP_7 PORTD7

// Base SPI @s for Input stage, channel and bias.
#define SPI_ADDR_INPUT_STAGE 0x00
#define SPI_ADDR_CHANNEL 0x08
#define SPI_ADDR_BIAS 0x20

// Output channel type:
#define VoLG 0
#define VoHG 1
#define VoSE 2
#define VoNC 3  //Not connected
#define VoANA 10  //MUSIC Output configured as analog
#define VoDIG 20  //MUSIC Output configured as digital

// ADC conversion factor 0~1023 -> 0 --> 3.3 V
#define V_LSB 0.0032258

// Target VDCs in LSBs:
#define VDC_TARGET_SE 0.3
#define MAX_VDC_DEVIATION 0.05
#define VDC_TARGET_SUM_P 0.9
#define VDC_TARGET_SUM_N 1.9
#define VCM_TARGET (VDC_TARGET_SUM_P + VDC_TARGET_SUM_N)/2
#define MAX_VDC_DEVIATION_SUM 0.01
#define VDC_MIN_VOLT_SWING_LSB 100

//EEPROM
#define EEPROM_CONFIG_START 65

//(See PCB schematic)
//MUSIC pin#    MUSIC pin name    ATmega pin#      ATmega pin name
//   48              VoSE0            28                ADC5
//   47              VoSE1            27                ADC4
//   46              VoSE2            25                ADC2
//   45              VoSE3            26                ADC3
//   36              VoSE4            23                ADC0
//   35              VoSE5            24                ADC1
//   34              VoSE6            22                ADC7
//   33              VoSE7            19                ADC6
//   43              VoHG+            26                ADC3
//   42              VoHG-            25                ADC2
//   39              VoLG+            24                ADC1
//   38              VoLG-            23                ADC0

const uint8_t chan_map_se[8] = { A5 , A4 , A2 , A3 , A0 , A1 , A7 , A6 };
const uint8_t chan_map_hg[2] = { A3 , A2 };
const uint8_t chan_map_lg[2] = { A1 , A0 };

const uint8_t resistor_map_se[8] = { 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 };
const uint8_t resistor_map_hg[2] = { 1 , 2 };
const uint8_t resistor_map_lg[2] = { 3 , 5 };

///////////////// BOARD OPERATION FUNCTIONS /////////////////

void MiniMusicBoard::init( bool debug_mode , uint8_t softspi_clk_div )
{
  Serial.print(F("INFO:  - Debug Mode: "));
  Serial.println(debug_mode ? "Yes" : "No" );

  this->debug_mode = debug_mode;
  mySPI.begin();
  mySPI.setClockDivider( softspi_clk_div );
  mySPI.setBitOrder(MSBFIRST);
  mySPI.setDataMode(SPI_MODE1);
  Serial.println(F("INFO: Initialization done!"));
}

void MiniMusicBoard::set_pin_modes( ){
  pinMode(ASIC_RST,OUTPUT);
}

void MiniMusicBoard::asic_reset(){
  digitalWrite( ASIC_RST , LOW );
  digitalWrite( ASIC_RST , HIGH );
  digitalWrite( ASIC_RST , LOW );
}

//Power Down FOREVER (unless reset or interrupt occurs):
void MiniMusicBoard::power_down( ){
  Serial.println(F("POWEROFF"));
  LowPower.powerDown(SLEEP_FOREVER, ADC_OFF, BOD_OFF);
}

/////////////////////////////////////////////////////////////

//////////////// USER INTERACTION FUNCTIONS /////////////////

bool MiniMusicBoard::wait_for_user_interaction( uint16_t milliseconds ){
  uint16_t i_milli = 0;
  bool data_available;
  Serial.println( F("WAITINGFORCOMMAND") );
  do{
    i_milli += 1;
    delay(1);
    data_available = (Serial.available()>1);
  }while( !data_available && i_milli<milliseconds );

  if( !data_available ){
    Serial.println( F("INFO: User interaction Timeout.") );
  }
  return data_available;
}

//Converts 2 ASCII characters into a byte value:
uint8_t MiniMusicBoard::get_serial_byte(){
  while( Serial.available()<2 );
  int high = toupper( Serial.read() );
  high = ( (high >= 'A') ? high - 'A' + 10 : high - '0' );
  int low  = toupper( Serial.read() );
  low = ( (low >= 'A') ? low - 'A' + 10 : low - '0' );
  return( ((uint8_t)high<<4) | ((uint8_t)low) );
}

/////////////////////////////////////////////////////////////

////////////////// EEPROM ACCESS FUNCTIONS //////////////////

void MiniMusicBoard::write_music_config_eeprom(){
  uint16_t addr, data_high, data_low;
  //First byte contains the number of SPI registers:
  uint8_t n_spi_registers = this->get_serial_byte();
  Serial.print(F("INFO: Writing "));
  Serial.print( n_spi_registers );
  Serial.println(F(" SPI Registers into the MiniMUSIC EEPROM:"));
  EEPROM.update( EEPROM_CONFIG_START , n_spi_registers );
  uint16_t i_eeprom_addr = EEPROM_CONFIG_START + 1;
  for( int i_reg=0 ; i_reg<n_spi_registers ; i_reg++ ){
    //3 Bytes: 1 for SPI Address, 2 for Data.
    addr      = this->get_serial_byte();
    data_high = this->get_serial_byte();
    data_low  = this->get_serial_byte();
    this->print_eeprom_spi_reg( i_eeprom_addr , addr , (data_high<<8)|data_low );
    this->eeprom_write_if_change( i_eeprom_addr , addr , data_high , data_low );
    i_eeprom_addr += 3;
    Serial.println(F("STORED"));
  }
  Serial.println(F("INFO: Data was written."));
}

void MiniMusicBoard::eeprom_write_if_change( uint16_t eeprom_addr , uint8_t spi_addr , uint8_t data_high , uint8_t data_low ){
  EEPROM.update( eeprom_addr   , spi_addr );
  EEPROM.update( eeprom_addr+1   , data_high );
  EEPROM.update( eeprom_addr+2   , data_low );
}

void MiniMusicBoard::read_music_config_eeprom( uint8_t music_pin, uint8_t *spi_addrs , uint16_t *spi_data )
{
  uint8_t addr, data_high, data_low;
  uint16_t i_eeprom_addr;
  if (music_pin == 0) i_eeprom_addr = 0x00;   // Bottom music chip
  if (music_pin == 1) i_eeprom_addr = 0x45;   // Top music chip

  for( int i_reg=0 ; i_reg < 23 ; i_reg++ ){
    addr      = EEPROM.read( i_eeprom_addr   );
    data_high = EEPROM.read( i_eeprom_addr+1 );
    data_low  = EEPROM.read( i_eeprom_addr+2 );
    i_eeprom_addr += 3;
    spi_addrs[i_reg] = addr;
    spi_data[i_reg] = (data_high<<8)|data_low;
  }
}

uint8_t MiniMusicBoard::read_music_serial(){
  uint8_t serial_len = EEPROM.read(0);
  Serial.print(F("INFO: Reading "));
  Serial.print( serial_len );
  Serial.println(F(" bytes serial"));
  if( serial_len == 0xFF ){
    Serial.println(F("ERROR: MiniMUSIC EEPROM is empty!"));
    Serial.println(F("SERIAL: 00"));
    this->erase_eeprom();
    return 0;
  }
  Serial.print(F("SERIAL: "));
  this->print_hex_2digits(serial_len);
  for( int i=0 ; i<serial_len; i++ ){
    this->print_hex_2digits( EEPROM.read( i + 1 ) );
  }
  Serial.println(F(""));
  return serial_len;
}

void MiniMusicBoard::write_music_serial(){
  //First byte contains the number of SPI registers:
  uint8_t serial_len = this->get_serial_byte();
  if( serial_len > EEPROM_CONFIG_START - 1 ) {
    serial_len = EEPROM_CONFIG_START -1 ;
  }
  Serial.print(F("INFO: Writing "));
  Serial.print( serial_len );
  Serial.println(F(" bytes serial into the MiniMUSIC EEPROM:"));
  EEPROM.write( 0 , serial_len );
  uint8_t d;
  for( int i=0 ; i<serial_len; i++ ){
    //3 Bytes: 1 for SPI Address, 2 for Data.
    d = this->get_serial_byte();
    EEPROM.write( i+1, d );
    this->print_hex_2digits(d);
    Serial.println(F(""));
    Serial.println(F("STORED"));
  }
  Serial.println(F("INFO: Data was written."));
}

void MiniMusicBoard::print_hex_2digits(int num) {
     char tmp[3];
     sprintf(tmp, "%02X", num);
     Serial.print(tmp);
}

void MiniMusicBoard::erase_eeprom(){
  Serial.println(F("INFO: Erasing EEPROM..."));
  int eeprom_length = EEPROM.length() - EEPROM_CONFIG_START;
  for (int i = EEPROM_CONFIG_START ; i < eeprom_length ; i++) {
    EEPROM.write(i, 0);
    if( i%64 == 0 ){
      Serial.print(F("INFO: "));
      Serial.print(i,DEC);
      Serial.print(F(" / "));
      Serial.println(eeprom_length,DEC);
    }
  }
  Serial.println(F("INFO: Done!"));
}

void MiniMusicBoard::print_eeprom_spi_reg( uint16_t eeprom_addr , uint16_t spi_addr , uint16_t spi_data ){
    //SPI Address:
    Serial.print(F("INFO: EEPROM[0x") );
    Serial.print( eeprom_addr , HEX );
    Serial.print(F(":0x") );
    Serial.print( eeprom_addr , HEX );
    Serial.print(F("]: SPI_ADDR = 0x") );
    Serial.println( spi_addr , HEX );
    //SPI Data:
    Serial.print(F("INFO: EEPROM[0x") );
    Serial.print( eeprom_addr+1 , HEX );
    Serial.print(F(":0x") );
    Serial.print( eeprom_addr+2 , HEX );
    Serial.print(F("]: SPI_DATA = 0x") );
    Serial.println( spi_data , HEX );
}

/////////////////////////////////////////////////////////////

/////////////////// SPI MASTER FUNCTIONS ////////////////////

void MiniMusicBoard::set_softspi_idle( uint8_t music_pin , bool idle ){
	//Slave select pin control:
	digitalWrite( music_pin, idle ? HIGH : LOW );
}

void MiniMusicBoard::spi_reg_write( uint8_t music_pin , uint8_t addr , uint16_t data , bool print_msg ){
  uint8_t data_msb, data_lsb, rx_msb, rx_lsb;
  uint16_t rx_data;
  data_msb = (uint8_t)(data>>8);
  data_lsb = (uint8_t)data;

  this->spi_write( music_pin , addr , data_msb , data_lsb );
  rx_data = this->spi_read( music_pin , addr , 0x0F , 0x0F );
  bool tx_ok = ( data == rx_data );
  if( print_msg ){
  	Serial.print(F("INFO: SPI @: 0x"));
  	Serial.print(addr,HEX);
  	Serial.print(F(" - Data: 0x"));
  	Serial.print(data,HEX);
    Serial.println( tx_ok ? " (OK)" : " (FAIL)" );
  }
  if( !tx_ok ){
    Serial.print(F("ERROR: Rx: 0x"));
    Serial.println(rx_data,HEX);
  }
}

uint16_t MiniMusicBoard::spi_read( uint8_t music_pin , uint8_t addr , uint8_t data_msb , uint8_t data_lsb ){
  digitalWrite(music_pin, LOW);
  mySPI.transfer( addr|0x80 );
   _delay_us(100);
  uint8_t rx_msb = mySPI.transfer( data_msb );
  _delay_us(100);
  uint8_t rx_lsb = mySPI.transfer( data_lsb );
  _delay_us(100);
  digitalWrite(music_pin, HIGH);
  uint16_t ret_val = (((uint16_t)rx_msb)<<8) | ((uint16_t)rx_lsb);
  if( this->debug_mode ) {
    Serial.print(F("DEBUG: Read SPI@0x"));
    Serial.print(addr,HEX);
    Serial.print(F(": 0x"));
    Serial.println(ret_val,HEX);
  }
  return ret_val;
}

void MiniMusicBoard::spi_write( uint8_t music_pin , uint8_t addr , uint8_t data_msb , uint8_t data_lsb ){
  uint16_t tx_val;
  //THIS SHOULD NOT OCCUR, BUT MODE1 DOES NOT WRITE DATA CORRECTLY!!!!
  //MUSICR1, MUSICR1a decodes address with the falling edge (as expected),
  //but stores data into the registers with the rising edge (bug).
  //To overcome this issue, SPI mode must be switched to Mode 0 during the data cycles.
  //mySPI.setDataMode(SPI_MODE1);
  digitalWrite(music_pin, LOW);
  mySPI.transfer( addr );
  _delay_us(100);
  mySPI.setDataMode(SPI_MODE0);
  mySPI.transfer( data_msb );
  _delay_us(100);
  mySPI.transfer( data_lsb );
  _delay_us(100);
  mySPI.setDataMode(SPI_MODE1);
  digitalWrite(music_pin, HIGH);
  if( this->debug_mode ) {
    tx_val = (((uint16_t)data_msb)<<8) | ((uint16_t)data_lsb);
    Serial.print(F("DEBUG: Write SPI@0x"));
    Serial.print(addr,HEX);
    Serial.print(F(": 0x"));
    Serial.println( tx_val );
  }
}

/////////////////////////////////////////////////////////////

/////////////// MUSIC THRESHOLD SCAN FUNCTIONS /////////////////

void MiniMusicBoard::threshold_scan( uint8_t music_pin , uint32_t n_meas_per_point ){
  uint16_t v_th_opt[8], eeprom_addr, eeprom_data;
  Serial.println(F("INFO: Starting Threshold Scan..."));
  //Gets the list of enabled channels (and digital). Test will only run on these channels:
  this->check_enabled_channels( music_pin , true );
  //Obtains the optimal VBG (all the enabled channels have a transition within this value):
  uint8_t opt_vbg = this->find_opt_vbg( music_pin , n_meas_per_point );
  Serial.print(F("RESULT: VBG = "));
  Serial.println(opt_vbg);
  //Per-channel threshold scan:
  for( int i_ch=0 ; i_ch<8 ; i_ch++ ){
    v_th_opt[i_ch] = 0;
    if( this->enabled_channel[i_ch] ){
      v_th_opt[i_ch] = this->threshold_scan_channel( music_pin , i_ch , n_meas_per_point );
      Serial.print(F("RESULT: V_TH["));
      Serial.print(i_ch);
      Serial.print(F("] = "));
      Serial.println(v_th_opt[i_ch]);
    }
  }
  Serial.println(F("RESULT: Command line options friendly copy&paste:"));
  Serial.print(F("RESULT: -g "));
  Serial.print(opt_vbg);
  Serial.print(F(" -v"));
  for( int i_ch=0 ; i_ch<8 ; i_ch++ ){
    Serial.print(F(" "));
    Serial.print(v_th_opt[i_ch]);
  }
  Serial.println("");
/*  Serial.println(F("INFO: Writing EEPROM With optimal VBG and Thresholds..."));
  //VBG
  eeprom_addr = 1+(SPI_ADDR_BIAS+2)*3;
  // BYTE_MSB: x x x D D D x x - BYTE_LSB: x x x x x x x x
  eeprom_data = (EEPROM.read( eeprom_addr+1 ) & 0x1C) | (opt_vbg<<2);
  EEPROM.write( eeprom_addr+1 , eeprom_data );
*/
}

uint8_t MiniMusicBoard::find_opt_vbg( uint8_t music_pin , uint32_t n_meas_per_point ){
  uint8_t i_vbg;
  uint16_t conf_mask;
  uint32_t count_ones_vth_000, count_ones_vth_511;
  bool found_vbg=false, transition_ok;
  //Initial VBG values:
  uint16_t conf_reg_vbg = this->spi_read( music_pin , SPI_ADDR_BIAS+2 , 0 , 0 );
  //Scans all VBG values trying to find the one which produces a transition in all MUSIC channels
  //when sweeping VTH from 0 to 511.
  Serial.println(F("INFO: Bandgap Voltage (VBG) scan..."));
  for( i_vbg=0 ; ((i_vbg<8)&&!found_vbg) ; i_vbg++ ){
    //Configuration of VBG:
    conf_reg_vbg = (conf_reg_vbg & 0xE3FF) | (((uint16_t)i_vbg)<<10);
    this->spi_reg_write( music_pin , SPI_ADDR_BIAS+2 , conf_reg_vbg , false );
    transition_ok = true;
    Serial.print(F("INFO: VBG = "));
    Serial.println(i_vbg);
    //For each channel we check VTH corners and checks if tarnsition occured:
    for( uint8_t i_ch=0 ; i_ch<8 ; i_ch++ ){
      if( this->enabled_channel[i_ch] ){
        conf_mask = this->spi_read( music_pin , SPI_ADDR_CHANNEL+i_ch , 0 , 0 );
        count_ones_vth_000 = this->get_counts_vth( music_pin , i_ch , conf_mask , 0   , n_meas_per_point );
        count_ones_vth_511 = this->get_counts_vth( music_pin , i_ch , conf_mask , 511 , n_meas_per_point );
        Serial.print(F("INFO: VoSE["));
        Serial.print(i_ch);
        Serial.print(F("]: V_TH=0: "));
        Serial.print((float)count_ones_vth_000/(float)n_meas_per_point);
        Serial.print(F(" --> V_TH=511: "));
        Serial.println((float)count_ones_vth_511/(float)n_meas_per_point);
        //If any of them is not fullfilled, transition_ok is false.
        transition_ok &= ( ( count_ones_vth_000<(n_meas_per_point>>1) ) && ( count_ones_vth_511>(n_meas_per_point>>1) ) );
      }
    }
    //When all transitions found, found_vbg goes 'true' and exits.
    found_vbg |= transition_ok;
  }
  Serial.print(F("INFO: Optimal Bandgap Voltage (VBG) is "));
  Serial.println(i_vbg-1);
  return (i_vbg-1);
}

uint16_t MiniMusicBoard::threshold_scan_channel( uint8_t music_pin , uint8_t channel , uint32_t n_meas_per_point ){
  uint32_t count_ones, count_ones_ant=0;
  uint16_t v_th_opt = 0;
  bool found=false;
  Serial.print( F("INFO: Running Threshold Scan on VoSE[") );
  Serial.print( channel , DEC );
  Serial.println( F("]...") );
  //Initial register value:
  uint16_t conf_mask = this->spi_read( music_pin , SPI_ADDR_CHANNEL+channel , 0 , 0 );
  //Iterates over the 512 thresholds, unless 50% duty cycle transition has been found.
  for( uint16_t i_th=0 ; (i_th<512)&&(count_ones!=n_meas_per_point) ; i_th++ ){
    count_ones = this->get_counts_vth( music_pin , channel , conf_mask , i_th , n_meas_per_point );
    if( (count_ones > (n_meas_per_point>>1)) && !found){
      v_th_opt = i_th;
      found = true;
    }
    //Serial.println( count_ones );
    if( (count_ones > 0) && (count_ones_ant == 0)){
      Serial.print( F("INFO: [000 .. ") );
      Serial.print( i_th-1 );
      Serial.println( F("] = 0.000") );
    }
    if( (count_ones > 0) && (count_ones < n_meas_per_point)){
      Serial.print( F("INFO: [") );
      Serial.print( i_th );
      Serial.print( F("] = ") );
      Serial.println( (float)count_ones/(float)n_meas_per_point , 3 );
    }
    if( (count_ones == n_meas_per_point) && (count_ones_ant < n_meas_per_point) ){
      Serial.print( F("INFO: [") );
      Serial.print( i_th );
      Serial.println( F(" .. 511] = 1.000") );
    }
    count_ones_ant = count_ones;
  }
  //Resets Vth = 0
  this->spi_reg_write( music_pin , SPI_ADDR_CHANNEL+channel , conf_mask & 0xF007 , false );
  Serial.print( F("INFO: Optimal V_TH for VoSE[") );
  Serial.print( channel , DEC );
  Serial.print( F("]: ") );
  Serial.println( v_th_opt , DEC );
  return v_th_opt;
}

//Configures a V_TH value to a given channel and counts the number of high levels of the digital signal.
uint32_t MiniMusicBoard::get_counts_vth( uint8_t music_pin , uint8_t channel , uint16_t conf_mask , uint16_t v_th , uint32_t n_meas_per_point ){
  conf_mask = (conf_mask & 0xF007) | (v_th<<3);
  this->spi_reg_write( music_pin , SPI_ADDR_CHANNEL+channel , conf_mask , false );
  uint32_t count_ones = 0;

  //Channels 6 & 7 are analog only! Compare reading with Vfs/2
  if( channel > 5 ){
    for( uint32_t i_meas=0 ; i_meas<n_meas_per_point ; i_meas++ ){
      count_ones += (uint32_t) (analogRead( chan_map_se[channel] ) > 511);
    }
  } else {
    for( uint32_t i_meas=0 ; i_meas<n_meas_per_point ; i_meas++ ){
      count_ones += (uint32_t) digitalRead( chan_map_se[channel] );
    }
  }
  return count_ones;
}

//Gets the list of enabled channels. Test will only run on these channels:
void MiniMusicBoard::check_enabled_channels( uint8_t music_pin , bool digital_ch_only  ){
  bool output_conf_ok;
  uint8_t ch_type;
  this->n_enabled_ch=0;
  Serial.println(F("RESULT: VoSE channel list:"));
  for( uint8_t i_ch=0 ; i_ch<8 ; i_ch++ ){
    //Checks that output is properly configurated (analog or digital).
    ch_type = this->check_output_type( music_pin , i_ch );
    output_conf_ok = (digital_ch_only && ch_type==VoDIG) || (!digital_ch_only && ch_type==VoANA);
    //Enable channel switch is the bit #10 of the input stage registers:
    this->enabled_channel[i_ch] = ( (this->spi_read( music_pin , SPI_ADDR_INPUT_STAGE+i_ch,0,0)&0x0400)==0x0400 ) && output_conf_ok;
    if( this->enabled_channel[i_ch] ){
      this->n_enabled_ch++;
    }
    Serial.print(F("RESULT:  - VoSE["));
    Serial.print(i_ch);
    Serial.print(F("]: "));
    if( ch_type==VoDIG ) Serial.print(F("Digital"));
    if( ch_type==VoANA ) Serial.print(F("Analog"));
    if( ch_type==VoNC  ) Serial.print(F("(Wrong Output Config)"));
    this->enabled_channel[i_ch] ? Serial.println(F(" (OK)")) : Serial.println(F(""));
  }
}

uint8_t MiniMusicBoard::check_output_type( uint8_t music_pin , uint8_t ch ){
  //Digital: when 'enDrvSE' (b12) = 0 and 'enCompSW' (b1) = 1
  //Analog: when 'enDrvSE' (b12) = 1
  uint16_t reg_val = this->spi_read( music_pin , SPI_ADDR_CHANNEL+ch,0,0);
  if( (reg_val&0x1002) == 0x0002){
    return VoDIG;
  } else if( (reg_val&0x1000) == 0x1000){
    return VoANA;
  } else {
    return VoNC;
  }
}


/////////////////////////////////////////////////////////////

/////////////// MUSIC VDC CALIBRATION FUNCTIONS /////////////////

//Configures VDC values: LG, HG, SE
void MiniMusicBoard::set_vdc( uint8_t music_pin , uint8_t value , uint8_t output_type ){
  uint8_t addr;
  if( output_type == VoSE ){
    addr = SPI_ADDR_BIAS+1;
  } else {
    addr = SPI_ADDR_BIAS;
  }
  //Reads initial configuration:
  uint16_t conf_mask = this->spi_read( music_pin , addr , 0 , 0 );
  if( output_type == VoHG ){
    conf_mask = (conf_mask & 0xFF00) | (uint16_t)value;
  } else {
    conf_mask = (conf_mask & 0x00FF) | ((uint16_t)value) << 8;
  }
  this->spi_reg_write( music_pin , addr , conf_mask , false );
}

//Configures VCM value
void MiniMusicBoard::set_vcm( uint8_t music_pin , uint8_t value ){
  uint8_t addr = SPI_ADDR_BIAS+1;
  //Reads initial configuration:
  uint16_t conf_mask = this->spi_read( music_pin , addr , 0 , 0 );
  conf_mask = (conf_mask & 0xFF00) | (uint16_t)value;
  this->spi_reg_write( music_pin , addr , conf_mask , false );
}

//This function checks the connectivity between MUSIC outputs (VoSE, VoHG, VoLG) and ATmega ADCs.
//Reminder: this uController has only 8 ADC channels, while MUSIC has 12 outputs. Signals are selected
//with an on-board 0-Ohm resistor.
void MiniMusicBoard::discover_vdc_outputs( uint8_t music_pin ){
  for( int i_ch=0 ; i_ch<8 ; i_ch++ ){
    //If channel is enabled and ADC is sensitive to changes in MUSIC VDC,
    //it means that the ADC channel is connected to a MUSIC SE output.
    if( this->enabled_channel[i_ch] ){
      Serial.print(F("INFO: Connectivity check with VoSE["));
      Serial.print(i_ch);
      Serial.println(F("]..."));
      //Even if the channel is enabled, we remove it from the channels to be tested.
      if( this->check_vdc_corners( music_pin , VoSE , chan_map_se[i_ch] ) == false ){
        this->enabled_channel[i_ch] = false;
        this->n_enabled_ch--;
      }
    } else {
      Serial.print(F("INFO: Connectivity check with VoSE["));
      Serial.print(i_ch);
      Serial.println(F("]: Disabled"));
    }
  }
  Serial.print(F("INFO: Number of VoSE connected channels: "));
  Serial.println(this->n_enabled_ch);
  Serial.println(F("INFO: Connectivity check with VoLG+..."));
  this->is_lg_enabled  = this->check_vdc_corners( music_pin , VoLG , chan_map_lg[0] );
  Serial.println(F("INFO: Connectivity check with VoLG-..."));
  this->is_lg_enabled &= this->check_vdc_corners( music_pin , VoLG , chan_map_lg[1] );
  Serial.println(F("INFO: Connectivity check with VoHG+..."));
  this->is_hg_enabled  = this->check_vdc_corners( music_pin , VoHG , chan_map_hg[0] );
  Serial.println(F("INFO: Connectivity check with VoHG-..."));
  this->is_hg_enabled &= this->check_vdc_corners( music_pin , VoHG , chan_map_hg[1] );
}

void MiniMusicBoard::print_connectivity_map(){
  Serial.println(F("RESULT: MiniMUSIC VDC connectivity map:"));
  for( int i_ch=0 ; i_ch<8 ; i_ch++ ){
    if( this->enabled_channel[i_ch] ){
      Serial.print(F("RESULT:  - VoSE["));
      Serial.print(i_ch);
      Serial.print(F("]: "));
      this->print_pin_correspondence( chan_map_se[i_ch] , resistor_map_se[i_ch] );
    }
  }
  if( this->is_lg_enabled ){
    Serial.print(F("RESULT:  - VoLG+  : "));
    this->print_pin_correspondence( chan_map_lg[0] , resistor_map_lg[0] );
    Serial.print(F("RESULT:  - VoLG-  : "));
    this->print_pin_correspondence( chan_map_lg[1] , resistor_map_lg[1] );
  }
  if( this->is_hg_enabled ){
    Serial.print(F("RESULT:  - VoHG+  : "));
    this->print_pin_correspondence( chan_map_hg[0] , resistor_map_hg[0] );
    Serial.print(F("RESULT:  - VoHG-  : "));
    this->print_pin_correspondence( chan_map_hg[1] , resistor_map_hg[1] );
  }
}

void MiniMusicBoard::print_pin_correspondence( uint8_t adc_ch_num , uint8_t resistor_num ){
    Serial.print(F("ADC["));
    Serial.print(adc_ch_num-A0);
    Serial.print(F("] (0-Ohm resistor on R"));
    Serial.print(resistor_num);
    Serial.println(F(")"));
}

bool MiniMusicBoard::check_vdc_corners( uint8_t music_pin , uint8_t output_type , uint8_t adc_pin_number ){
  uint16_t vo_vdc_min, vo_vdc_max;
  bool is_connected;
  this->set_vdc( music_pin , 0xFF , output_type );     //Highest VDC = 0xFF
  delay(10);                               //MUSIC outputs take some time to change.
  vo_vdc_min = analogRead( adc_pin_number );
  this->set_vdc( music_pin , 0x00 , output_type );     //Lowest VDC = 0xFF
  delay(10);
  vo_vdc_max = analogRead( adc_pin_number );
  Serial.print(F("INFO: Vo(min) = "));
  Serial.print((float)vo_vdc_min*V_LSB);
  Serial.println(F(" V"));
  Serial.print(F("INFO: Vo(max) = "));
  Serial.print((float)vo_vdc_max*V_LSB);
  Serial.println(F(" V"));
  //VDC_MIN_VOLT_SWING_LSB -> > x Volts of difference between max and min VDC.
  //If condition is fullfilled, the ADC is connected to a MUSIC output.
  is_connected = ( vo_vdc_max - vo_vdc_min > VDC_MIN_VOLT_SWING_LSB );
  is_connected ? Serial.println(F("INFO: Connected")) : Serial.println(F("INFO: Disconnected"));
  //Insufficient voltage swing, but > 0: something wrong may happen
  if( !is_connected && (vo_vdc_max - vo_vdc_min > VDC_MIN_VOLT_SWING_LSB/2) ){
    Serial.print(F("WARNING: Non-zero voltage swing, but lower than "));
    Serial.print((float)VDC_MIN_VOLT_SWING_LSB*V_LSB);
    Serial.println(F(" V"));
  }
  return is_connected;
}

//Read N times the analog voltage correspondence.
float MiniMusicBoard::vdc_meas_avg( uint8_t adc_pin_number , uint16_t n_avg_meas ){
  uint32_t meas_acc_lsb = 0;
  float meas_value;
  for( uint16_t i_meas=0 ; i_meas<n_avg_meas ; i_meas++ ){
    meas_acc_lsb += (uint32_t)analogRead( adc_pin_number );
  }
  meas_value = ( (float)meas_acc_lsb * V_LSB / (float)n_avg_meas );
  return meas_value;
}

uint8_t MiniMusicBoard::vdc_se_calib( uint8_t music_pin , uint16_t n_avg_meas ){
  uint8_t i_vdc=0;
  bool err_sign=true, err_sign_ant=true;
  float vdc_meas_volts, vdc_max, vdc_min, vdc_acc, err_v=1000, err_ant_v;
  //Loops until the measured value is in the error range or lowest error found.
  do{
    err_sign_ant = err_sign;
    err_ant_v = err_v;
    vdc_max = 0;
    vdc_min = 10;
    vdc_acc = 0;
    this->set_vdc( music_pin , i_vdc , VoSE );  //Changes in VDC takes some time to be measurable
    delay(1);
    for( int i_ch=0 ; i_ch<8 ; i_ch++ )
    {
      //For each connected channel, VDC is measured, accumulates and checks if it's a maximum or minimum.
      if( this->enabled_channel[i_ch] ){
        vdc_meas_volts = this->vdc_meas_avg( chan_map_se[i_ch] , n_avg_meas );
        if( vdc_max < vdc_meas_volts ){
          vdc_max = vdc_meas_volts;
        }
        if( vdc_min > vdc_meas_volts ){
          vdc_min = vdc_meas_volts;
        }
        vdc_acc += vdc_meas_volts;
      }
    }
    //Computes the average error.
    err_v = (vdc_acc/(float)this->n_enabled_ch) - VDC_TARGET_SE;
    err_sign = (err_v > 0);
    Serial.print(F("INFO: VDC_CH = "));
    Serial.print(i_vdc);
    Serial.print(F(" LSB, VoSE (min,max,avg) = "));
    Serial.print(vdc_min,3);
    Serial.print(F(" , "));
    Serial.print(vdc_max,3);
    Serial.print(F(" , "));
    Serial.print((vdc_acc/(float)this->n_enabled_ch),4);
    Serial.print(F(" V. Error = "));
    Serial.print(round(err_v*1000));
    Serial.print(F(" mV (Target = "));
    Serial.print(VDC_TARGET_SE);
    Serial.println(F(" V)."));
    err_v = abs( err_v );
    //To speedup the procedure:
    if( err_v > 1 ){
      i_vdc += 25;
    } else if (err_v > 0.3 ) {
      i_vdc += 5;
    } else {
      i_vdc++;
    }
  }while( err_sign == err_sign_ant );
  //Returns the optimal VDC.
  return ( i_vdc-2 );
}

// With the maximum differential voltage swing (bias regs = 0xFF) VCM register
// is swept until the measured VCM is in the middle (VCM_TARGET = 1.4 volts).
uint8_t MiniMusicBoard::vcm_calib( uint8_t music_pin , uint16_t n_avg_meas ){
  float vdc_lg[2], vdc_hg[2], meas_vcm, err_v=1000, err_ant_v;
  uint8_t i_vcm=0;
  this->set_vdc( music_pin , 0xFF , VoLG );
  this->set_vdc( music_pin , 0xFF , VoHG );
  do{
    err_ant_v = err_v;
    this->set_vcm( music_pin , i_vcm );  //Changes in VDC takes some time to be measurable
    delay(1);
    vdc_lg[0] = this->vdc_meas_avg( chan_map_lg[0] , n_avg_meas );
    vdc_lg[1] = this->vdc_meas_avg( chan_map_lg[1] , n_avg_meas );
    vdc_hg[0] = this->vdc_meas_avg( chan_map_hg[0] , n_avg_meas );
    vdc_hg[1] = this->vdc_meas_avg( chan_map_hg[1] , n_avg_meas );
    if( this->is_lg_enabled && this->is_hg_enabled ){
      meas_vcm = (vdc_lg[0] + vdc_lg[1] + vdc_hg[0] + vdc_hg[1]) / 4.0;
    } else if( this->is_hg_enabled ){
      meas_vcm = (vdc_hg[0] + vdc_hg[1]) / 2.0;
    } else if( this->is_hg_enabled ){
      meas_vcm = (vdc_lg[0] + vdc_lg[1]) / 2.0;
    }
    err_v = abs( meas_vcm-VCM_TARGET );
    Serial.print(F("INFO: VCM = "));
    Serial.print(i_vcm);
    Serial.print(F(" LSB, meas = "));
    Serial.print(meas_vcm,3);
    Serial.print(F(" V. Error = "));
    Serial.print(round(err_v*1000));
    Serial.print(F(" mV (Target = "));
    Serial.print(VCM_TARGET);
    Serial.println(F(" V)."));
    if( this->debug_mode ){
      Serial.print(F("DEBUG: LG+ (V) = "));
      Serial.println(vdc_lg[0],3);
      Serial.print(F("DEBUG: LG- (V) = "));
      Serial.println(vdc_lg[1],3);
      Serial.print(F("DEBUG: HG+ (V) = "));
      Serial.println(vdc_hg[0],3);
      Serial.print(F("DEBUG: HG- (V) = "));
      Serial.println(vdc_hg[1],3);
    }
    //To speedup the procedure:
    if( err_v > 0.2 ){
      i_vcm += 10;
    } else if (err_v > 0.1 ) {
      i_vcm += 5;
    } else {
      i_vcm++;
    }
  }while( err_v<err_ant_v );
  return ( i_vcm-2 );
}

uint8_t MiniMusicBoard::vdc_sum_calib( uint8_t music_pin , uint8_t output_type, uint16_t n_avg_meas ){
  float vdc_sum[2],  err_v=1000, err_ant_v;
  uint8_t i_vdc_sum=0;
  do{
    err_ant_v = err_v;
    this->set_vdc( music_pin , i_vdc_sum , output_type );  //Changes in VDC takes some time to be measurable
    delay(1);
    if( output_type == VoLG ){
      Serial.print(F("INFO: VDC_LG = "));
      vdc_sum[0] = this->vdc_meas_avg( chan_map_lg[0] , n_avg_meas );
      vdc_sum[1] = this->vdc_meas_avg( chan_map_lg[1] , n_avg_meas );
    } else {
      Serial.print(F("INFO: VDC_HG = "));
      vdc_sum[0] = this->vdc_meas_avg( chan_map_hg[0] , n_avg_meas );
      vdc_sum[1] = this->vdc_meas_avg( chan_map_hg[1] , n_avg_meas );
    }
    err_v = max( abs(vdc_sum[0]-VDC_TARGET_SUM_P) , abs(vdc_sum[1]-VDC_TARGET_SUM_N) );
    Serial.print(i_vdc_sum);
    Serial.print(F(" LSB. Error = "));
    Serial.print(round(err_v*1000));
    Serial.print(F(" mV (Targets = "));
    Serial.print(VDC_TARGET_SUM_P);
    Serial.print(F(", "));
    Serial.print(VDC_TARGET_SUM_N);
    Serial.println(F(" V)."));
    if( this->debug_mode ){
      Serial.print(F("DEBUG: SUM+ (V) = "));
      Serial.println(vdc_sum[0],3);
      Serial.print(F("DEBUG: SUM- (V) = "));
      Serial.println(vdc_sum[1],3);
    }
    //To speedup the procedure:
    if( err_v > 0.5 ){
      i_vdc_sum += 50;
    } else if (err_v > 0.1 ) {
      i_vdc_sum += 5;
    } else {
      i_vdc_sum++;
    }
  }while( err_v<err_ant_v );
  return ( i_vdc_sum-2 );
}

void MiniMusicBoard::vdc_calib( uint8_t music_pin , uint16_t n_avg_meas, bool store_eeprom ){
  uint8_t vdc_se_opt, vcm_opt, vdc_lg_opt, vdc_hg_opt;
  //fills 'is_se_ch_connected' boolean array and counts how many enabled chanels are. Only channels configured as analog.
  this->check_enabled_channels( music_pin , false );
  //Gets MUSIC --> ATmega VDC outputs connectivity.
  this->discover_vdc_outputs( music_pin );
  this->print_connectivity_map();
  if( this->n_enabled_ch > 0 ){
    //Calibrates only the discovered VoSE channels:
    Serial.println(F("INFO: Calibrating MUSIC Single-Ended VDC Outputs..."));
    vdc_se_opt = this->vdc_se_calib( music_pin , n_avg_meas );
    Serial.print(F("INFO: Optimal VDC_CH = "));
    Serial.println(vdc_se_opt);
    Serial.print(F("RESULT: VDC_CH = "));
    Serial.println(vdc_se_opt);
    if(store_eeprom){this->set_vdc_eeprom( vdc_se_opt , VoSE );}
    this->set_vdc( music_pin , vdc_se_opt , VoSE );
  }
  if( this->is_lg_enabled || this->is_hg_enabled ){
    //Calibrates VCM. Target VCM = (1.9 + 0.9)/2 = 1.4 Volts:
    Serial.println(F("INFO: Calibrating MUSIC Differential sum common mode (VCM)..."));
    vcm_opt = this->vcm_calib( music_pin , n_avg_meas );
    Serial.print(F("INFO: Optimal VCM = "));
    Serial.println(vcm_opt);
    Serial.print(F("RESULT: VCM = "));
    Serial.println(vcm_opt);
    if(store_eeprom){this->set_vcm_eeprom( vcm_opt );}    //EEPROM storage
    this->set_vcm( music_pin , vcm_opt );           //SPI MUSIC register change.
  }
  if( this->is_lg_enabled ){
    //Calibrates VDC_LG: finds the value which minimizes the error.
    //Expected VoLG+ = 1.9 Volts, VoLG- = 0.9 Volts
    Serial.println(F("INFO: Calibrating MUSIC Low-Gain Differential sum (VoLG)..."));
    vdc_lg_opt = this->vdc_sum_calib( music_pin , VoLG , n_avg_meas );
    Serial.print(F("INFO: Optimal VDC_LG = "));
    Serial.println(vdc_lg_opt);
    Serial.print(F("RESULT: VDC_LG = "));
    Serial.println(vdc_lg_opt);
    if(store_eeprom){this->set_vdc_eeprom( vdc_lg_opt , VoLG );}
    this->set_vdc( music_pin , vdc_lg_opt , VoLG );
  }
  if( this->is_hg_enabled ){
    Serial.println(F("INFO: Calibrating MUSIC Low-Gain Differential sum (VoHG)..."));
    vdc_hg_opt = this->vdc_sum_calib( music_pin , VoHG , n_avg_meas );
    Serial.print(F("INFO: Optimal VDC_HG = "));
    Serial.println(vdc_hg_opt);
    Serial.print(F("RESULT: VDC_HG = "));
    Serial.println(vdc_hg_opt);
    if(store_eeprom){this->set_vdc_eeprom( vdc_hg_opt , VoHG );}
    this->set_vdc( music_pin , vdc_hg_opt , VoHG );
  }
}

//Configures VDC values: LG, HG, SE
void MiniMusicBoard::set_vdc_eeprom( uint8_t value , uint8_t output_type ){
  Serial.println("INFO: Updating EEPROM with optimal VDC...");
  uint8_t eeprom_addr = EEPROM_CONFIG_START;
  switch( output_type ){
    case VoSE: eeprom_addr += 0x35;  break;
    case VoLG: eeprom_addr += 0x32;  break;
    case VoHG: eeprom_addr += 0x33;  break;
  }
  Serial.print(F("INFO: EEPROM[0x") );
  Serial.print( eeprom_addr , HEX );
  Serial.print(F("] = 0x") );
  Serial.println( value , HEX );
  EEPROM.update( eeprom_addr , value );
}

//Configures VCM value
void MiniMusicBoard::set_vcm_eeprom( uint8_t value ){
  Serial.println("INFO: Updating EEPROM with optimal VCM...");
  uint8_t eeprom_addr = 0x36 + EEPROM_CONFIG_START;
  Serial.print(F("INFO: EEPROM[0x") );
  Serial.print( eeprom_addr , HEX );
  Serial.print(F("] = 0x") );
  Serial.println( value , HEX );
  EEPROM.update( eeprom_addr , value );
}

void MiniMusicBoard::write_to_spi( uint8_t music_pin ) {
  uint16_t addr, data_high, data_low;
  //First byte contains the number of SPI registers:
  uint8_t n_spi_registers = this->get_serial_byte();
  Serial.print(F("INFO: Writing "));
  Serial.print( n_spi_registers );
  Serial.println(F(" SPI Registers into the MUSIC"));
  set_softspi_idle( music_pin , true );
  for( int i_reg=0 ; i_reg<n_spi_registers ; i_reg++ ){
    //3 Bytes: 1 for SPI Address, 2 for Data.
    addr      = this->get_serial_byte();
    data_high = this->get_serial_byte();
    data_low  = this->get_serial_byte();
    spi_reg_write(music_pin , addr, (data_high<<8)|data_low );
    Serial.println(F("STORED"));
  }
  Serial.println(F("INFO: Data was written."));
}

void MiniMusicBoard::read_from_spi( uint8_t music_pin ) {
  uint8_t d1, d2;
  uint16_t data;
  Serial.print(F("AUTO: 17"));
  for( uint8_t i = 0; i < 8; i++ ) {
    data = spi_read(music_pin,i,d1,d2);
    this->print_hex_2digits( i );
    this->print_hex_2digits( data>>8 );
    this->print_hex_2digits( data&0xFF );
    data = spi_read(music_pin,i+8,d1,d2);
    this->print_hex_2digits( i+8 );
    this->print_hex_2digits( data>>8 );
    this->print_hex_2digits( data&0xFF );
  }
  for( uint8_t i = 0x20; i < 0x27; i++ ) {
    data = spi_read(music_pin,i,d1,d2);
    this->print_hex_2digits( i );
    this->print_hex_2digits( data>>8 );
    this->print_hex_2digits( data&0xFF );
  }
  Serial.println(F(""));
}

void MiniMusicBoard::configure_music( uint8_t music_pin, uint8_t spi_addrs[], uint16_t spi_data[])
{
  for( int i_reg=0 ; i_reg < 23 ; i_reg++ )
  {    
    spi_reg_write( music_pin , spi_addrs[i_reg] , spi_data[i_reg] , true);
  }
}
