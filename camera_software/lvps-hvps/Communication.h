#ifndef _COMMUNICATION_H_
#define _COMMUNICATION_H_

#include "SPI.h"

unsigned char SPI_Init();
void SPI_Read(unsigned char* data);
void SPI_Write(unsigned char data[3]);

#endif
