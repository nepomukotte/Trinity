/* =============================================================================
 * ControllerDeviceIsp.hpp
 */

#ifndef get_zboard_hw_ControllerDeviceIsp_hpp_INCLUDED
#define get_zboard_hw_ControllerDeviceIsp_hpp_INCLUDED

#include "mdaq/hw/AbstractControllerDevice.hpp"
#include "mdaq/hw/HwControllerException.h"
#include <fcntl.h>
#include <sys/mman.h>
#include <map>
#include <iostream>
#include <cstdlib>

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
class ControllerDeviceIsp : public AbstractControllerDevice<T, P>
{
public:
	typedef std::map<P, T> MemoryMap;

	ControllerDeviceIsp(uint64_t& baseAddress);
	virtual ~ControllerDeviceIsp();

	void declareRegister(P address);
	void undeclareRegister(P address);
	T read(P address);
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
ControllerDeviceIsp<T, P>::ControllerDeviceIsp(uint64_t& baseAddress) 
: deviceAddress(baseAddress), dip_mapped_size(0), AbstractControllerDevice<T,P>(baseAddress) 
{
	dip_memfd = open("/dev/mem", O_RDWR);
	if (dip_memfd == -1)
	{
	  throw("ControllerDeviceIsp : Can't open /dev/mem");
	  std::cout << "ControllerDeviceIsp : Can't open /dev/mem" << std::endl;
	}
	else
	{
	  std::cout << "ControllerDeviceIsp 0x" << std::hex << baseAddress << " : /dev/mem opened" << std::endl;
	//*(P( (char *)(baseAddress) + 0x4 )) = 0x00000010;  // REG1 = 0x10
  	//*(P( (char *)(baseAddress) )) = 0x10000000;
	
	}

}

/**
* destructor
*/
template<typename T, typename P>
ControllerDeviceIsp<T, P>::~ControllerDeviceIsp() 
{
	if (munmap(dip_mapped_base, dip_mapped_size) == -1) {
	  throw("ControllerDeviceIsp : Can't unmap memory from user space");
	  std::cout << "ControllerDeviceIsp : Can't unmap memory from user space" << std::endl;
	}
	close(dip_memfd);
	std::cout << "ControllerDeviceIsp 0x" << deviceAddress << " : /dev/mem is closed" << std::endl;
}

//____________________________________________________________________________//

/**
 * Declares a register at the specified address.
 */
template<typename T, typename P>
void ControllerDeviceIsp<T, P>::declareRegister(P address)
{
	// Map one page of memory into user space such that the device is in that page, but it may not
	// be at the start of the page

	if ( ((int)(address) + sizeof(T)) > dip_mapped_size )
	{
		uint64_t mask = (uint64_t)address + sizeof(T) - 1;
		int newSize   = (int)(address) + sizeof(T);
		if (dip_mapped_size != 0)
			if (munmap(dip_mapped_base, dip_mapped_size) == -1) {
			  throw("ControllerDeviceIsp : Can't unmap memory from user space");
			  std::cout << "ControllerDeviceIsp : Can't unmap memory from user space" << std::endl;
			}
		dip_mapped_base = mmap(0, newSize, PROT_READ | PROT_WRITE, MAP_SHARED, dip_memfd, deviceAddress & ~mask);
		dip_mapped_size = newSize;
		if (dip_mapped_base == (void *) -1) 
		{
		  throw("ControllerDeviceIsp : Can't map the memory to user space.");
		  std::cout << "ControllerDeviceIsp : Can't map the memory to user space." << std::endl;
		}
		else
		{
		  //std::cout << "ControllerDeviceIsp : Memory mapped at address : " << dip_mapped_base << ", with size = " << dip_mapped_size << std::endl; 
		}
	}
}

/**
 * Undeclare the register at the specified address.
 * Further reading and writing of this address will result to an exception.
 */
template<typename T, typename P>
void ControllerDeviceIsp<T, P>::undeclareRegister(P address)
{
}



/**
 * Returns the value of the register at the specified address.
 */
template<typename T, typename P>
T ControllerDeviceIsp<T, P>::read(P address)
{
	
 	return *(P( (char*)(dip_mapped_base) + (int)address ));
  
}

/**
 * Sets the value of the register at the specified address.
 */
template<typename T, typename P>
void ControllerDeviceIsp<T, P>::write(P address, T value)
{
	*(P( (char *)(dip_mapped_base) + (int)address )) = (int)value;
 	
}

} // namespace hw
} // namespace zboard
} // namespace get

#endif // get_zboard_hw_ControllerDeviceIsp_hpp_INCLUDED
