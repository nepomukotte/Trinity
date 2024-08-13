/* =============================================================================
 * ControllerDeviceAdc.hpp
 */

#ifndef get_zboard_hw_ControllerDeviceAdc_hpp_INCLUDED
#define get_zboard_hw_ControllerDeviceAdc_hpp_INCLUDED

#include "mdaq/hw/AbstractControllerDevice.hpp"
#include "mdaq/hw/HwControllerException.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <map>
#include <iostream>
#include <cstdlib>

#define REG0 0
#define REG1 4
#define REG2 8
#define REG3 12
#define REG4 16
#define REG5 20
#define REG6 24
#define REG7 28
#define REG8 32
#define REG9 36

using mdaq::hw::AbstractControllerDevice;

namespace get
{
namespace zboard
{
namespace hw
{
/**
 * This class implements a device policy that simulates reading and writing
 * registers into an addressed space.
 */
template<typename T, typename P = T*>
class ControllerDeviceAdc : public AbstractControllerDevice<T, P>
{
public:
	typedef std::map<P, T> MemoryMap;

	ControllerDeviceAdc(uint64_t& baseAddress);
	virtual ~ControllerDeviceAdc();

	void declareRegister(P address);
	void undeclareRegister(P address);
	T read(P address);
	T get_address(P address);
	T get_reg_read(P address);
	T get_reg_write(P address);
	void write(P address, T value);
private:
	MemoryMap mem;
	uint64_t deviceAddress;
	int dip_memfd;
	void *dip_mapped_base;
	int dip_mapped_size;
};

/**
* constructor
*/
template<typename T, typename P>
ControllerDeviceAdc<T, P>::ControllerDeviceAdc(uint64_t& baseAddress) 
: deviceAddress(baseAddress), dip_mapped_size(0), AbstractControllerDevice<T,P>(baseAddress) 
{
	dip_memfd = open("/dev/mem", O_RDWR);
	if (dip_memfd == -1)
	{
	  throw("ControllerDeviceAdc : Can't open /dev/mem");
	  std::cout << "ControllerDeviceAdc : Can't open /dev/mem" << std::endl;
	}
	else
	{
	  std::cout << "ControllerDeviceAdc 0x" << std::hex << baseAddress << " : /dev/mem opened" << std::endl;
	//*(P( (char *)(baseAddress) + 0x4 )) = 0x00000010;  // REG1 = 0x10
  	//*(P( (char *)(baseAddress) )) = 0x10000000;
	
	}

}

/**
* destructor
*/
template<typename T, typename P>
ControllerDeviceAdc<T, P>::~ControllerDeviceAdc() 
{
	if (munmap(dip_mapped_base, dip_mapped_size) == -1) {
	  throw("ControllerDeviceAdc : Can't unmap memory from user space");
	  std::cout << "ControllerDeviceAdc : Can't unmap memory from user space" << std::endl;
	}
	close(dip_memfd);
	std::cout << "ControllerDeviceAdc 0x" << deviceAddress << " : /dev/mem is closed" << std::endl;
}

//____________________________________________________________________________//

/**
 * Declares a register at the specified address.
 */
template<typename T, typename P>
void ControllerDeviceAdc<T, P>::declareRegister(P address)
{
	// Map one page of memory into user space such that the device is in that page, but it may not
	// be at the start of the page

	if ( ((int)(address) + sizeof(T)) > dip_mapped_size )
	{
		uint64_t mask = (uint64_t)address + sizeof(T) - 1;
		int newSize   = (int)(address) + sizeof(T);
		if (dip_mapped_size != 0)
			if (munmap(dip_mapped_base, dip_mapped_size) == -1) {
			  throw("ControllerDeviceAdc : Can't unmap memory from user space");
			  std::cout << "ControllerDeviceAdc : Can't unmap memory from user space" << std::endl;
			}
		dip_mapped_base = mmap(0, newSize, PROT_READ | PROT_WRITE, MAP_SHARED, dip_memfd, deviceAddress & ~mask);
		dip_mapped_size = newSize;
		if (dip_mapped_base == (void *) -1) 
		{
		  throw("ControllerDeviceAdc : Can't map the memory to user space.");
		  std::cout << "ControllerDeviceAdc : Can't map the memory to user space." << std::endl;
		}
		else
		{
		  //std::cout << "ControllerDeviceAdc : Memory mapped at address : " << dip_mapped_base << ", with size = " << dip_mapped_size << std::endl; 
		}
	}
}

/**
 * Undeclare the register at the specified address.
 * Further reading and writing of this address will result to an exception.
 */
template<typename T, typename P>
void ControllerDeviceAdc<T, P>::undeclareRegister(P address)
{
}

/**
 * Decode  REGX  for write.
 */
template<typename T, typename P>
T ControllerDeviceAdc<T, P>::get_reg_write(P address)
{
int regX = 0;

switch((int)address) 
{
case 0:
	regX = REG0;
	break;
case 1:
	regX = REG1;
	break;
case 2:
	regX = REG2;
	break;
case 3:
	regX = REG3;
	break;
case 4:
	regX = REG4;
	break;
case 5:
	regX = REG5;
	break;
case 6:
	regX = REG6;
	break;
case 7:
	regX = REG7;
	break;
case 8:
	regX = REG8;
	break;
case 9:
	regX = REG9;
	break;					
}
return regX ;
}
/**
 * Returns the value of the register at the specified address.
 */
template<typename T, typename P>
T ControllerDeviceAdc<T, P>::read(P address)
{

 return *(P( (char*)(dip_mapped_base) + (int)address ));
  
}

/**
 * Sets the value of the register at the specified address.
 */
template<typename T, typename P>
void ControllerDeviceAdc<T, P>::write(P address, T value)
{
  int reg = get_reg_write(address); 
  *(P( (char *)(dip_mapped_base) + reg )) = (int)value;
  int val = (int)address << 4 ;
  val = val +  0x80000000 ;
  *(P( (char *)(dip_mapped_base) + REG9)) = val;
  
   std::cout << "ADC access : write at adress: " << address << " with value: " << value << std::endl;
}

} // namespace hw
} // namespace zboard
} // namespace get

#endif // get_zboard_hw_ControllerDeviceAdc_hpp_INCLUDED
