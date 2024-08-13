#ifndef _MUSIC_CHIP_H_
#define _MUSIC_CHIP_H_

#include <Arduino.h>
#include "SoftSPI.h"

#define music_top_rst   A0    // PC0
#define music_bot_rst   A1    // PC1


class MiniMusicBoard
{
  bool debug_mode = false;
  bool enabled_channel[8];
  uint8_t n_enabled_ch;
  bool is_lg_enabled , is_hg_enabled;

  void     print_eeprom_spi_reg( uint16_t eeprom_addr , uint16_t spi_addr , uint16_t spi_data );
  void     print_hex_2digits(int num);
  SoftSPI mySPI;
    
  public:

  MiniMusicBoard() : mySPI(musicChipMOSI, musicChipMISO, musicChipSCLK){};
  void init( bool debug_mode, uint8_t softspi_clk_div );
  void set_pin_modes();
  void asic_reset();
  bool wait_for_user_interaction( uint16_t milliseconds );
  uint8_t get_serial_byte();
  void write_music_config_eeprom();
  void read_music_config_eeprom( uint8_t music_pin, uint8_t *spi_addrs , uint16_t *spi_data );
  void eeprom_write_if_change( uint16_t eeprom_addr , uint8_t spi_addr , uint8_t data_high , uint8_t data_low );
  void erase_eeprom();
  void power_down( );

  void set_softspi_idle( uint8_t music_pin , bool idle );
  void spi_reg_write( uint8_t music_pin , uint8_t addr , uint16_t data , bool print_msg=true );

  void threshold_scan( uint8_t music_pin , uint32_t n_meas_per_point );
  uint8_t find_opt_vbg( uint8_t music_pin , uint32_t n_meas_per_point );
  uint16_t threshold_scan_channel( uint8_t music_pin , uint8_t channel , uint32_t n_meas_per_point );
  uint32_t get_counts_vth( uint8_t music_pin , uint8_t channel , uint16_t conf_mask , uint16_t v_th , uint32_t n_meas_per_point );
  void check_enabled_channels( uint8_t music_pin , bool digital_ch_only );
  uint8_t check_output_type( uint8_t music_pin , uint8_t ch );
  void set_vdc( uint8_t music_pin , uint8_t value , uint8_t output_type );
  void set_vcm( uint8_t music_pin , uint8_t value );
  void discover_vdc_outputs( uint8_t music_pin );
  void print_connectivity_map();
  void print_pin_correspondence( uint8_t adc_ch_num , uint8_t resistor_num  );
  bool check_vdc_corners( uint8_t music_pin , uint8_t output_type , uint8_t adc_pin_number );
  float vdc_meas_avg( uint8_t adc_pin_number , uint16_t n_avg_meas );
  uint8_t vdc_se_calib( uint8_t music_pin , uint16_t n_avg_meas );
  uint8_t vcm_calib( uint8_t music_pin , uint16_t n_avg_meas );
  uint8_t vdc_sum_calib( uint8_t music_pin , uint8_t output_type, uint16_t n_avg_meas );
  void vdc_calib( uint8_t music_pin , uint16_t n_avg_meas, bool store_eeprom );
  void set_vdc_eeprom( uint8_t value , uint8_t output_type );
  void set_vcm_eeprom( uint8_t value );  

  uint8_t read_music_serial();
  void write_music_serial();

  void     spi_write( uint8_t music_pin , uint8_t addr , uint8_t data_msb , uint8_t data_lsb );
  uint16_t spi_read ( uint8_t music_pin , uint8_t addr , uint8_t data_msb , uint8_t data_lsb );

  void write_to_spi( uint8_t music_pin );
  void read_from_spi( uint8_t music_pin );
  void configure_music (uint8_t music_pin, uint8_t spi_addrs[], uint16_t spi_data[]);
};

#endif
