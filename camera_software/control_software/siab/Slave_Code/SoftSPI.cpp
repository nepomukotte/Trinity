
#include "SoftSPI.h"

SoftSPI::SoftSPI(uint8_t musicChipMOSI, uint8_t musicChipMISO, uint8_t musicChipSCLK) {
    _musicChipMOSI = musicChipMOSI;
    _musicChipMISO = musicChipMISO;
    _musicChipSCLK = musicChipSCLK;
    _musicdelay = 2;
    _musiccke = 0;
    _musicckp = 0;
    _musicorder = MSBFIRST;
}

void SoftSPI::begin() {
    pinMode(_musicChipMOSI, OUTPUT);
    pinMode(_musicChipMISO, INPUT);
    pinMode(_musicChipSCLK, OUTPUT);
}

void SoftSPI::end() {
    pinMode(_musicChipMOSI, INPUT);
    pinMode(_musicChipMISO, INPUT);
    pinMode(_musicChipSCLK, INPUT);
}

void SoftSPI::setBitOrder(uint8_t musicorder) {
    _musicorder = musicorder & 1;
}

void SoftSPI::setDataMode(uint8_t mode) {
    switch (mode) {
        case SPI_MODE0:
            _musicckp = 0;
            _musiccke = 0;
            break;
        case SPI_MODE1:
            _musicckp = 0;
            _musiccke = 1;
            break;
        case SPI_MODE2:
            _musicckp = 1;
            _musiccke = 0;
            break;
        case SPI_MODE3:
            _musicckp = 1;
            _musiccke = 1;
            break;
    }
}

void SoftSPI::setClockDivider(uint8_t div) {
    switch (div) {
        case SPI_CLOCK_DIV2:
            _musicdelay = 2;
            break;
        case SPI_CLOCK_DIV4:
            _musicdelay = 4;
            break;
        case SPI_CLOCK_DIV8:
            _musicdelay = 8;
            break;
        case SPI_CLOCK_DIV16:
            _musicdelay = 16;
            break;
        case SPI_CLOCK_DIV32:
            _musicdelay = 32;
            break;
        case SPI_CLOCK_DIV64:
            _musicdelay = 64;
            break;
        case SPI_CLOCK_DIV128:
            _musicdelay = 128;
            break;
        default:
            _musicdelay = 128;
            break;
    }
}

uint8_t SoftSPI::transfer(uint8_t val) {
    uint8_t out = 0;
    if (_musicorder == MSBFIRST) {
        uint8_t v2 = 
            ((val & 0x01) << 7) |
            ((val & 0x02) << 5) |
            ((val & 0x04) << 3) |
            ((val & 0x08) << 1) |
            ((val & 0x10) >> 1) |
            ((val & 0x20) >> 3) |
            ((val & 0x40) >> 5) |
            ((val & 0x80) >> 7);
        val = v2;
    }

    uint8_t del = _musicdelay >> 1;

    uint8_t bval = 0;
    //if CPHA=1, flip clock line

    if(_musiccke) // if CPHA=1 , this sets the clock before the first sampling
    {
        digitalWrite(_musicChipSCLK, _musicckp ? LOW : HIGH);
    }
    
    for (uint8_t bit = 0; bit < 8; bit++) {

        //setup the bit on MOSI
        digitalWrite(_musicChipMOSI, val & (1<<bit) ? HIGH : LOW);
        
        //delay 1
        for (uint8_t i = 0; i < del; i++) {
            asm volatile("nop");
        }

        // sample here  
        //read MISO, write MOSI
        bval = digitalRead(_musicChipMISO);
        if (_musicorder == MSBFIRST) {
            out <<= 1;
            out |= bval;
        } else {
            out >>= 1;
            out |= bval << 7;
        }
       
        // flip clock, but this should be the sampling edge
        if(_musiccke)
        {
          digitalWrite(_musicChipSCLK, _musicckp ? HIGH : LOW);
        } else {
          digitalWrite(_musicChipSCLK, _musicckp ? LOW : HIGH);
        }
       
        //delay 2 
        for (uint8_t i = 0; i < del; i++) {
            asm volatile("nop");
        }

        // flip clock, but this should be the non-sampling edge
        if(!(_musiccke && bit==7))
        {
            if(_musiccke)
            {
              digitalWrite(_musicChipSCLK, _musicckp ? LOW : HIGH);
            } else {
              digitalWrite(_musicChipSCLK, _musicckp ? HIGH : LOW);
            }
        } 
    }
    
    return out;
}

uint16_t SoftSPI::transfer16(uint16_t data)
{
	union {
		uint16_t val;
		struct {
			uint8_t lsb;
			uint8_t msb;
		};
	} in, out;
  
	in.val = data;

	if ( _musicorder == MSBFIRST ) {
		out.msb = transfer(in.msb);
		out.lsb = transfer(in.lsb);
	} else {
		out.lsb = transfer(in.lsb);
		out.msb = transfer(in.msb);
	}

	return out.val;
}
