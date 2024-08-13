/* =============================================================================
 * AGetControllerDevice.hpp
 * -----------------------------------------------------------------------------
 * template class mdaq::hw::AGetControllerDevice
 * Created on: 23 September 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of
 * the Mordicus Real-Time Software Development Framework.
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#ifndef get_drv_AGetControllerDevice_hpp_INCLUDED
#define get_drv_AGetControllerDevice_hpp_INCLUDED

#include "get/drv/AGetControllerDevice.h"
#include "mdaq/utl/CmdException.h"
#include "utl/Logging.h"
#include <map>
#include <iomanip>

#define THROW_CMD_EXCEPTION(msg)  { \
	LOG_ERROR() << msg;                  \
	throw ::mdaq::utl::CmdException(msg);\
}

namespace get
{
namespace drv
{

template<typename T, typename P>
AGetControllerDevice< T, P >::AGetControllerDevice(uint64_t & baseAddress)
	: mdaq::hw::AbstractControllerDevice<T,P>(baseAddress),
	  deviceAddress(reinterpret_cast<char*>(baseAddress))
{
	lengthPtr_uint32        = reinterpret_cast<Length_uint32*>    (deviceAddress + lengthOffset);
	readWritePtr_uint32     = reinterpret_cast<ReadWrite_uint32*> (deviceAddress + readWriteOffset);
	regAddressPtr_uint32    = reinterpret_cast<RegAddress_uint32*> (deviceAddress + regAddressOffset);
	statusPtr_uint32        = reinterpret_cast<Status_uint32*>    (deviceAddress + statusOffset);
	dataInPtr_uint64        = reinterpret_cast<BodyWrite_uint64*> (deviceAddress + dataInOffset);
	dataOutPtr_uint64       = reinterpret_cast<BodyRead_uint64*>  (deviceAddress + dataOutOffset);
	/* initialization */
	regAddressPtr_uint32->setValue(0);
	lengthPtr_uint32->setValue(0);
	dataInPtr_uint64[B_127_64].setValue(0);
	dataInPtr_uint64[B_63_0].setValue(0);
	dataOutPtr_uint64[B_127_64].setValue(0);
	dataOutPtr_uint64[B_63_0].setValue(0);
	fictiveAddress_uint32 = 0;
	/* association between the number of the register and its length */
	ASSOCIATE_LENGTH_BIT(reg0);
	ASSOCIATE_LENGTH_BIT(reg1);
	ASSOCIATE_LENGTH_BIT(reg2);
	ASSOCIATE_LENGTH_BIT(reg3);
	ASSOCIATE_LENGTH_BIT(reg4);
	ASSOCIATE_LENGTH_BIT(reg5);
	ASSOCIATE_LENGTH_BIT(reg6);
	ASSOCIATE_LENGTH_BIT(reg7);
	ASSOCIATE_LENGTH_BIT(reg8); //lsb of reg8
	ASSOCIATE_LENGTH_BIT(reg9); //lsb of reg9
	ASSOCIATE_LENGTH_BIT(reg10);
	ASSOCIATE_LENGTH_BIT(reg11);
	ASSOCIATE_LENGTH_BIT(reg12);
	ASSOCIATE_LENGTH_BIT(reg8msb); //msb of reg8
	ASSOCIATE_LENGTH_BIT(reg9msb); //msb of reg9
	ASSOCIATE_LENGTH_BIT(hit);     //lsb of hit register
	ASSOCIATE_LENGTH_BIT(hitmsb);  //msb of  register
}


//____________________________________________________________________________//

/**
 * Declares a register at the specified address.
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::declareRegister(P /*address*/)
{
}

/**
 * Undeclare the register at the specified address.
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::undeclareRegister(P /*address*/)
{
}

/**
 * Returns the length of the data associated with the address
 */
template<typename T, typename P>
uint32_t AGetControllerDevice<T, P>::length(P address)
{
	typename MemoryMap::const_iterator itReg = mem.find(reinterpret_cast<uint32_t>(address));
	if(itReg == mem.end())
	{
		throw typename mdaq::hw::AbstractControllerDevice<T, P>::UndeclaredRegister(address);
	}
	return itReg->second;
}

/**
 * prints the value of the data register and of other meaningful registers
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::debug()
{
	uint32_t* dataIn32;
	uint32_t* dataOut32;

	dataIn32 = reinterpret_cast<uint32_t*>(dataInPtr_uint64);
	dataOut32 = reinterpret_cast<uint32_t*>(dataOutPtr_uint64);

	std::cout << std::hex << std::noshowbase << "======== AGET" << std::endl;
	std::cout << "regAddress    = 0x" << regAddressPtr_uint32->value() << std::endl;
	std::cout << "lengthAddress = " << std::dec << lengthPtr_uint32->value() << std::hex << std::endl;
	std::cout << "dataIn        = 0x" << std::flush;
	for (int i = 0; i < 4; i++)
	{
		std::cout << std::hex << std::setfill('0') << std::setw(8) << *(dataIn32+i);
	}
	std::cout << std::endl << "dataOut       = 0x" << std::flush;
	for (int i = 0; i < 4; i++)
	{
		std::cout << std::hex << std::setfill('0') << std::setw(8) << *(dataOut32+i);
	}
	std::cout << std::endl << std::dec << std::setfill(' ');
}

/**
 * Sends an exception if the register is write-only
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::checkIfReadingIsAllowed(P address)
{
	if (reinterpret_cast<uint32_t>(address) == reg0)
		THROW_CMD_EXCEPTION("AGet Register is Write-only");
}

/**
 * Sends an exception if the register is read-only
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::checkIfWritingIsAllowed(P address)
{
	if (reinterpret_cast<uint32_t>(address) == reg5)
		THROW_CMD_EXCEPTION("AGet Register is Read-only");
}

/**
 * Sends an exception if the register is read-only
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::setRegAddress(P address)
{
	fictiveAddress_uint32 = reinterpret_cast<uint32_t>(address);

	switch (fictiveAddress_uint32)
	{
	case reg8msb:
		regAddressPtr_uint32->setValue(reg8);
		break;
	case reg9msb:
		regAddressPtr_uint32->setValue(reg9);
		break;
	case reg0: case reg1: case reg2: case reg3: case reg4:
	case reg5: case reg6: case reg7: case reg8: case reg9:
	case reg10: case reg11: case reg12: case hit : case hitmsb :
		regAddressPtr_uint32->setValue(fictiveAddress_uint32);
		break;
	default:
		THROW_CMD_EXCEPTION("AGet Register does not exist");
		break;
	}
}

/**
* Writes on the length register the length of the register
*/
template<typename T, typename P>
void AGetControllerDevice<T, P>::setLength(P address)
{
	lengthPtr_uint32->setValue(length(address));
}

/**
 * Sends read command.
 *
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::sendReadCommand()
{
	assertNotBusy();
	readWritePtr_uint32->setValue(1);
}

/**
 * Sends write command.
 *
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::sendWriteCommand()
{
	assertNotBusy();
	readWritePtr_uint32->setValue(0);
}

/**
 * Waits until the AGET controller is not busy.
 * @param numTrials Maximum number of checks of the busy signal.
 * @exception mdaq::utl::CmdException AGET controller is busy.
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::assertNotBusy(size_t numTrials)
{
	while ( (statusPtr_uint32->busy() == 1) && (numTrials > 0))
	{
		if (statusPtr_uint32->error() == 1)
		{
			THROW_CMD_EXCEPTION("AGET controller error");
		}
		numTrials--;
	}
	if (numTrials == 0)
	{
		unsigned* ptr = (unsigned*) statusPtr_uint32;
		std::cout << "=>>> Value @ 0x" << std::hex << ptr << " = 0x" << ptr[0] << std::dec << std::endl;
		THROW_CMD_EXCEPTION("AGET controller is busy");
	}
}

/**
 * Waits for the busy signal of the state machine to go down
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::checkEndOfCommand()
{
	if (statusPtr_uint32->error() == 1)
	{
		THROW_CMD_EXCEPTION("AGET controller error");
	}
	else
	{
		assertNotBusy();
	}
}


/**
 * Fills the body part of the protocol for a reading.
 * Here, it is just a reinitialisation of dataIn and dataOut
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::processBodyForReading()
{
	dataInPtr_uint64[B_127_64].setValue(0);
	dataInPtr_uint64[B_63_0].setValue(0);
	dataOutPtr_uint64[B_127_64].setValue(0);
	dataOutPtr_uint64[B_63_0].setValue(0);
}

/**
 * Fills the body part of the protocol for a writing
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::processBodyForWriting(T value)
{

//	cout << "======== AGET (bodyForWriting)" << endl;

	dataInPtr_uint64[B_127_64].setValue(0);
	dataInPtr_uint64[B_63_0].setValue(0);
	dataOutPtr_uint64[B_127_64].setValue(0);
	dataOutPtr_uint64[B_63_0].setValue(0);

	switch (fictiveAddress_uint32)
	{
	case reg0:
		dataInPtr_uint64[B_127_64].set_reg0(value);
		break;
	case reg1:
		dataInPtr_uint64[B_127_64].set_reg1(value);
		break;
	case reg2:
		dataInPtr_uint64[B_127_64].set_reg2(value);
		break;
	case reg3:
		dataInPtr_uint64[B_127_64].set_reg3(value);
		break;
	case reg4:
		dataInPtr_uint64[B_127_64].set_reg4(value);
		break;
	case reg5:
		dataInPtr_uint64[B_127_64].set_reg5(value);
		break;
	case reg6:
		dataInPtr_uint64[B_127_64].set_reg6(value);
		break;
	case reg7:
		dataInPtr_uint64[B_127_64].set_reg7(value);
		break;
	case hit:
		read(reinterpret_cast<P>(hit));
		{
			BodyRead_uint64 value_uint64;
			value_uint64.setValue(value);
			dataInPtr_uint64[B_127_64].set_msb4(dataOutPtr_uint64[B_127_64].lsb4());
			dataInPtr_uint64[B_127_64].set_lsb60(value_uint64.msb60());
			dataInPtr_uint64[B_63_0].set_msb4(value_uint64.lsb4());
		}
		break;
	case hitmsb:
		read(reinterpret_cast<P>(hit));
		dataInPtr_uint64[B_127_64].set_msb4(value);
		dataInPtr_uint64[B_127_64].set_lsb60(dataOutPtr_uint64[B_63_0].msb60());
		dataInPtr_uint64[B_63_0].set_msb4(dataOutPtr_uint64[B_63_0].lsb4());
		break;
	case reg8:
		read(reinterpret_cast<P>(reg8));
		dataInPtr_uint64[B_127_64].setValue(dataOutPtr_uint64[B_127_64].value());
		dataInPtr_uint64[B_63_0].setValue(value);
		break;
	case reg8msb:
		read(reinterpret_cast<P>(reg8));
		dataInPtr_uint64[B_127_64].setValue(value);
		dataInPtr_uint64[B_63_0].setValue(dataOutPtr_uint64[B_63_0].value());
		break;
	case reg9:
		read(reinterpret_cast<P>(reg9));
		dataInPtr_uint64[B_127_64].setValue(dataOutPtr_uint64[B_127_64].value());
		dataInPtr_uint64[B_63_0].setValue(value);
		break;
	case reg9msb:
		read(reinterpret_cast<P>(reg9));
		dataInPtr_uint64[B_127_64].setValue(value);
		dataInPtr_uint64[B_63_0].setValue(dataOutPtr_uint64[B_63_0].value());
		break;
	case reg10:
		dataInPtr_uint64[B_127_64].set_reg10(value);
		break;
	case reg11:
		dataInPtr_uint64[B_127_64].set_reg11(value);
		break;
	case reg12:
		dataInPtr_uint64[B_127_64].set_reg12(value);
		break;
	default:
		THROW_CMD_EXCEPTION("AGet Register does not exist");
		break;
	}

}

/**
 * Reads the correct part of the 128-bit body
 */
template<typename T, typename P>
T AGetControllerDevice<T, P>::processBodyAfterReading()
{
	T value;

	switch (fictiveAddress_uint32)
	{
	case reg8msb: case reg9msb :
		value = dataOutPtr_uint64[B_127_64].value();
		break;
	case hitmsb :
		value = dataOutPtr_uint64[B_127_64].lsb4();
		break;
	default:
		value = dataOutPtr_uint64[B_63_0].value();
		break;
	}

	return value;
}

/**
 * Returns the value of the register at the specified address.
 */
template<typename T, typename P>
T AGetControllerDevice<T, P>::read(P address)
{
	T value;
	checkIfReadingIsAllowed(address);
	setRegAddress(address);
	setLength(address);
	processBodyForReading();
	sendReadCommand();
	checkEndOfCommand();
	value = processBodyAfterReading();
//	debug();

	return value;
}

/**
 * Fills the body part of the protocol for a writing
 */
template<typename T, typename P>
void AGetControllerDevice<T, P>::write(P address, T value)
{
	checkIfWritingIsAllowed(address);
	setRegAddress(address);
	setLength(address);
	processBodyForWriting(value);
	sendWriteCommand();
	checkEndOfCommand();
//	debug();
}

}
}

#undef B_127_64
#undef B_63_0
#undef ASSOCIATE_LENGTH_BIT

#endif // get_drv_AGetControllerDevice_hpp_INCLUDED
