#include "Communication.h"

unsigned char SPI_Init()
{
  SPI.begin();
  return 0;
}

void SPI_Read(unsigned char* data)
{
  SPI.beginTransaction(SPISettings(8000000,MSBFIRST,SPI_MODE2));
  digitalWrite(SS,LOW);
  data[0]=SPI.transfer(0x00);
  data[1]=SPI.transfer(0x00);
  data[2]=SPI.transfer(0x00);
  digitalWrite(SS,HIGH);
  SPI.endTransaction();
  SPI.beginTransaction(SPISettings(8000000,MSBFIRST,SPI_MODE3));
}

void SPI_Write(unsigned char data[3])
{
  SPI.beginTransaction(SPISettings(8000000,MSBFIRST,SPI_MODE0));
  digitalWrite(SS,LOW);
  SPI.transfer(data[0]);
  SPI.transfer(data[1]);
  SPI.transfer(data[2]);
  SPI.transfer(data[0]);
  SPI.transfer(data[1]);
  SPI.transfer(data[2]);
  digitalWrite(SS,HIGH);
  SPI.endTransaction();
  SPI.beginTransaction(SPISettings(8000000,MSBFIRST,SPI_MODE3));
}
