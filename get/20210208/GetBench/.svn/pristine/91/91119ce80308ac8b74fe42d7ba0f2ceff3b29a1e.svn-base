/* =============================================================================
 * ControllerDeviceAget.hpp
 */

#ifndef get_zboard_hw_ControllerDeviceAget_hpp_INCLUDED
#define get_zboard_hw_ControllerDeviceAget_hpp_INCLUDED

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
class ControllerDeviceAget : public AbstractControllerDevice<T, P>
{
public:
	typedef std::map<P, T> MemoryMap;

	ControllerDeviceAget(uint64_t& baseAddress);
	virtual ~ControllerDeviceAget();

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
ControllerDeviceAget<T, P>::ControllerDeviceAget(uint64_t& baseAddress) 
: deviceAddress(baseAddress), dip_mapped_size(0), AbstractControllerDevice<T,P>(baseAddress) 
{
	dip_memfd = open("/dev/mem", O_RDWR);
	if (dip_memfd == -1)
	{
	  throw("ControllerDeviceAget : Can't open /dev/mem");
	  std::cout << "ControllerDeviceAget : Can't open /dev/mem" << std::endl;
	}
	else
	{
	  std::cout << "ControllerDeviceAget 0x" << std::hex << baseAddress << " : /dev/mem opened" << std::endl;
	//*(P( (char *)(baseAddress) + 0x4 )) = 0x00000010;  // REG1 = 0x10
  	//*(P( (char *)(baseAddress) )) = 0x10000000;
	
	}

}

/**
* destructor
*/
template<typename T, typename P>
ControllerDeviceAget<T, P>::~ControllerDeviceAget() 
{
	if (munmap(dip_mapped_base, dip_mapped_size) == -1) {
	  throw("ControllerDeviceAget : Can't unmap memory from user space");
	  std::cout << "ControllerDeviceAget : Can't unmap memory from user space" << std::endl;
	}
	close(dip_memfd);
	std::cout << "ControllerDeviceAget 0x" << deviceAddress << " : /dev/mem is closed" << std::endl;
}

//____________________________________________________________________________//

/**
 * Declares a register at the specified address.
 */
template<typename T, typename P>
void ControllerDeviceAget<T, P>::declareRegister(P address)
{
	// Map one page of memory into user space such that the device is in that page, but it may not
	// be at the start of the page

	if ( ((int)(address) + sizeof(T)) > dip_mapped_size )
	{
		uint64_t mask = (uint64_t)address + sizeof(T) - 1;
		int newSize   = (int)(address) + sizeof(T);
		if (dip_mapped_size != 0)
			if (munmap(dip_mapped_base, dip_mapped_size) == -1) {
			  throw("ControllerDeviceAget : Can't unmap memory from user space");
			  std::cout << "ControllerDeviceAget : Can't unmap memory from user space" << std::endl;
			}
		dip_mapped_base = mmap(0, newSize, PROT_READ | PROT_WRITE, MAP_SHARED, dip_memfd, deviceAddress & ~mask);
		dip_mapped_size = newSize;
		if (dip_mapped_base == (void *) -1) 
		{
		  throw("ControllerDeviceAget : Can't map the memory to user space.");
		  std::cout << "ControllerDeviceAget : Can't map the memory to user space." << std::endl;
		}
		else
		{
		  //std::cout << "ControllerDeviceAget : Memory mapped at address : " << dip_mapped_base << ", with size = " << dip_mapped_size << std::endl; 
		}
	}
}

/**
 * Undeclare the register at the specified address.
 * Further reading and writing of this address will result to an exception.
 */
template<typename T, typename P>
void ControllerDeviceAget<T, P>::undeclareRegister(P address)
{
}
/**
 * Decode  the new address of register.
 */
template<typename T, typename P>
T ControllerDeviceAget<T, P>::get_address(P address)
{
int  new_adress ;

switch((int)address) 
{
case 0:
case 1:
case 2:
case 3:
case 4:
case 5:
case 6:
case 7:
case 8:
case 9:
case 10:
case 11:
case 12: 
	new_adress = (int)address ;	
	break;
case 13: 
	new_adress = 3;	
	break;
case 14: 
	new_adress = 4;	
	break;
case 15: 
	new_adress = 6;	
	break;
case 16: 
	new_adress = 7;	
	break;
		
case 17: 
	new_adress = 8;	
	break;
case 18: 
	new_adress = 9;	
	break;	
case 19: 
	new_adress = 10;	
	break;
case 20: 
	new_adress = 11;	
	break;
case 21: 
	new_adress = 8;	
	break;
case 22: 
	new_adress = 9;	
	break;
case 23: 
	new_adress = 8;	
	break;	
case 24: 
	new_adress = 9;	
	break;					
default :
	break;
					
}
return new_adress ;
}

/**
 * Decode  REGX  for read.
 */
template<typename T, typename P>
T ControllerDeviceAget<T, P>::get_reg_read(P address)
{
int regX = 0;

switch((int)address) 
{
case 0:
case 1:
case 2:
case 3:
case 4:
case 5:
case 6:
case 7:
case 8:
case 9:
case 10:
case 11:
case 12: 
	regX = REG9;
	break;
case 13: 
case 14: 
case 15: 
case 16: 		
case 17: 
case 18: 
case 19: 
case 20: 
	regX = REG8;
	break;
case 21: 
case 22: 
	regX = REG7;;
	break;
case 23: 
case 24: 
	regX = REG6;;
	break;					
default :
	break;
					
}
return regX ;
}
/**
 * Decode  REGX  for write.
 */
template<typename T, typename P>
T ControllerDeviceAget<T, P>::get_reg_write(P address)
{
int regX = 0;

switch((int)address) 
{
case 0:
case 1:
case 2:
case 3:
case 4:
case 5:
case 6:
case 7:
case 8:
case 9:
case 10:
case 11:
case 12: 
	regX = REG2;
	break;
case 13: 
case 14: 
case 15: 
case 16: 		
case 17: 
case 18: 
case 19: 
case 20: 
	regX = REG3;
	break;
case 21: 
case 22: 
	regX = REG4;;
	break;
case 23: 
case 24: 
	regX = REG5;;
	break;					
default :
	break;
					
}
return regX ;
}
/**
 * Returns the value of the register at the specified address.
 */
template<typename T, typename P>
T ControllerDeviceAget<T, P>::read(P address)
{
// décodage adresse
  int  new_adress = get_address(address);
  int reg = get_reg_read(address);
  //*(P( (char *)(dip_mapped_base) + 0x4 )) = (int)address + 0x00000180;
  *(P( (char *)(dip_mapped_base) + 0x4 )) = (int)new_adress + 0x00000180;
  *(P( (char *)(dip_mapped_base) )) = 0x10000000; // send tree pulses for slow mode
  std::cout << "AGET Slow access : read at adress: " << address << std::endl;
  *(P( (char *)(dip_mapped_base) )) = 0x80000000;
  int val = (int) *(P( (char*)(dip_mapped_base)));
  std::cout << "Frame sent ..." << std::endl;
  uint32_t lsb = *(P( (char*)(dip_mapped_base) + reg));
  return lsb ;
  //return *(P( (char*)(dip_mapped_base) + 0x24));
 // return *(P( (char*)(dip_mapped_base) + (int)address ));
  
}

/**
 * Sets the value of the register at the specified address.
 */
template<typename T, typename P>
void ControllerDeviceAget<T, P>::write(P address, T value)
{
  int  new_adress = get_address(address);
  int reg = get_reg_write(address);
 // *(P( (char *)(dip_mapped_base) + (int)address )) = value;
  *(P( (char *)(dip_mapped_base) + REG1 )) = (int)new_adress + 0x100;
  *(P( (char *)(dip_mapped_base) + reg )) = (int)value;
  *(P( (char *)(dip_mapped_base) )) = 0x10000000;
   std::cout << "AGET Slow access : write at adress: " << address << " with value: " << value << std::endl;
  *(P( (char *)(dip_mapped_base) )) = 0x80000000;
}

} // namespace hw
} // namespace zboard
} // namespace get

#endif // get_zboard_hw_ControllerDeviceAget_hpp_INCLUDED
