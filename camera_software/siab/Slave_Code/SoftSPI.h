
#ifndef _SOFTSPI_H
#define _SOFTSPI_H

#define DEBUG_MODE false
#define SOFTSPI_CLOCK_DIVIDER SPI_CLOCK_DIV2

#include <SPI.h>

const int musicChipSCLK = 5;      // PD5
const int musicChipMISO = 6;      // PD6
const int musicChipMOSI = 7;      // PD7
const int musicChipCSBot = 8;     // PB0
const int musicChipCSTop = 9;     // PB1


class SoftSPI : public SPIClass {
    private:
        uint8_t _musiccke;
        uint8_t _musicckp;
        uint8_t _musicdelay;
        uint8_t _musicChipMISO;
        uint8_t _musicChipMOSI;
        uint8_t _musicChipSCLK;
        uint8_t _musicorder;

    public:
        SoftSPI(uint8_t musicChipMOSI, uint8_t musicChipMISO, uint8_t musicChipSCLK);
        void begin();
        void end();
        void setBitOrder(uint8_t);
        void setDataMode(uint8_t);
        void setClockDivider(uint8_t);
        uint8_t transfer(uint8_t);
		    uint16_t transfer16(uint16_t data);
};

#endif
