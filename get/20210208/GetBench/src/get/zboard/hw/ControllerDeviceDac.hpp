/* =============================================================================
 * ControllerDeviceDac.hpp
 */

#ifndef get_zboard_hw_ControllerDeviceDac_hpp_INCLUDED
#define get_zboard_hw_ControllerDeviceDac_hpp_INCLUDED

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
class ControllerDeviceDac : public AbstractControllerDevice<T, P>
{
public:
	typedef std::map<P, T> MemoryMap;

	ControllerDeviceDac(uint64_t& baseAddress);
	virtual ~ControllerDeviceDac();

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
ControllerDeviceDac<T, P>::ControllerDeviceDac(uint64_t& baseAddress) 
: deviceAddress(baseAddress), dip_mapped_size(0), AbstractControllerDevice<T,P>(baseAddress) 
{
	dip_memfd = open("/dev/mem", O_RDWR);
	if (dip_memfd == -1)
	{
	  throw("ControllerDeviceDac : Can't open /dev/mem");
	  std::cout << "ControllerDeviceDac : Can't open /dev/mem" << std::endl;
	}
	else
	{
	  std::cout << "ControllerDeviceDac 0x" << std::hex << baseAddress << " : /dev/mem opened" << std::endl;
	//*(P( (char *)(baseAddress) + 0x4 )) = 0x00000010;  // REG1 = 0x10
  	//*(P( (char *)(baseAddress) )) = 0x10000000;
	
	}

}

/**
* destructor
*/
template<typename T, typename P>
ControllerDeviceDac<T, P>::~ControllerDeviceDac() 
{
	if (munmap(dip_mapped_base, dip_mapped_size) == -1) {
	  throw("ControllerDeviceDac : Can't unmap memory from user space");
	  std::cout << "ControllerDeviceDac : Can't unmap memory from user space" << std::endl;
	}
	close(dip_memfd);
	std::cout << "ControllerDeviceDac 0x" << deviceAddress << " : /dev/mem is closed" << std::endl;
}

//____________________________________________________________________________//

/**
 * Declares a register at the specified address.
 */
template<typename T, typename P>
void ControllerDeviceDac<T, P>::declareRegister(P address)
{
	// Map one page of memory into user space such that the device is in that page, but it may not
	// be at the start of the page

	if ( ((int)(address) + sizeof(T)) > dip_mapped_size )
	{
		uint64_t mask = (uint64_t)address + sizeof(T) - 1;
		int newSize   = (int)(address) + sizeof(T);
		if (dip_mapped_size != 0)
			if (munmap(dip_mapped_base, dip_mapped_size) == -1) {
			  throw("ControllerDeviceDac : Can't unmap memory from user space");
			  std::cout << "ControllerDeviceDac : Can't unmap memory from user space" << std::endl;
			}
		dip_mapped_base = mmap(0, newSize, PROT_READ | PROT_WRITE, MAP_SHARED, dip_memfd, deviceAddress & ~mask);
		dip_mapped_size = newSize;
		if (dip_mapped_base == (void *) -1) 
		{
		  throw("ControllerDeviceDac : Can't map the memory to user space.");
		  std::cout << "ControllerDeviceDac : Can't map the memory to user space." << std::endl;
		}
		else
		{
		  //std::cout << "ControllerDeviceDac : Memory mapped at address : " << dip_mapped_base << ", with size = " << dip_mapped_size << std::endl; 
		}
	}
}

/**
 * Undeclare the register at the specified address.
 * Further reading and writing of this address will result to an exception.
 */
template<typename T, typename P>
void ControllerDeviceDac<T, P>::undeclareRegister(P address)
{
}

/**
 * Decode  Address.
 */
template<typename T, typename P>
T ControllerDeviceDac<T, P>::get_address(P address)
{
int new_address = 0;

switch((int)address) 
{
case 0:
	new_address = 0x0000;
	break;
case 1:
	new_address = 0x0200;
	break;
case 2:
	new_address = 0x0D00;
	break;
case 3:
	new_address = 0x0E00;
	break;
case 4:
	new_address = 0x0F00;
	break;
case 5:
	new_address = 0x1000;
	break;
case 6:
	new_address = 0x1100;
	break;
case 7:
	new_address = 0x1400;
	break;			
}
return new_address;
}


/**
 * Returns the value of the register at the specified address.
 */
template<typename T, typename P>
T ControllerDeviceDac<T, P>::read(P address)
{
	int new_address = get_address(address);
	int val = 0x8000 + new_address ;
	
	*(P( (char *)(dip_mapped_base) + REG1 )) = (int)val;
	*(P( (char *)(dip_mapped_base) + REG0)) = 0x80000000;
	std::cout << "DAC access : read at adress: " << address  << std::endl;
 	return *(P( (char*)(dip_mapped_base) + REG2 ));
  
}

/**
 * Sets the value of the register at the specified address.
 */
template<typename T, typename P>
void ControllerDeviceDac<T, P>::write(P address, T value)
{
  	int new_address = get_address(address);
	int val = new_address + (int)value ;
	*(P( (char *)(dip_mapped_base) + REG1 )) = (int)val;
	*(P( (char *)(dip_mapped_base) + REG0)) = 0x80000000;
	std::cout << "DAC access : write at adress: " << address  << " with value: " << value << std::endl;
 	
}

} // namespace hw
} // namespace zboard
} // namespace get

#endif // get_zboard_hw_ControllerDeviceDac_hpp_INCLUDED
