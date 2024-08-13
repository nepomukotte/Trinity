/* =============================================================================
 * AsAdControllerDevice.hpp
 * -----------------------------------------------------------------------------
 * template class mdaq::hw::AsAdControllerDevice
 * Created on: 23 September 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Joël Chavas (joel.chavas@cea.fr)
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

#ifndef get_drv_AsAdControllerDevice_hpp_INCLUDED
#define get_drv_AsAdControllerDevice_hpp_INCLUDED

#include "AsAdControllerDevice.h"
#include <iostream>

namespace get
{
namespace drv
{

/*===============================================================
 * Template Class implementation
 *==============================================================*/

template<typename T, typename P>
AsAdControllerDevice<T,P>::AsAdControllerDevice(uint64_t& baseAddress)
	: mdaq::hw::AbstractControllerDevice<T,P>(baseAddress),
	  deviceAddress(reinterpret_cast<char*>(baseAddress))
{
	moduleAddressPtr_uint32 = reinterpret_cast<ModuleAddress_uint32*> (deviceAddress + moduleAddressOffset);
	lengthPtr_uint32        = reinterpret_cast<Length_uint32*>    (deviceAddress + lengthOffset);
	readWritePtr_uint32     = reinterpret_cast<ReadWrite_uint32*> (deviceAddress + readWriteOffset);
	dataInPtr_uint64        = reinterpret_cast<BodyWrite_uint64*> (deviceAddress + dataInOffset);
	dataOutPtr_uint64       = reinterpret_cast<BodyRead_uint64*>  (deviceAddress + dataOutOffset);
	statusPtr_uint32        = reinterpret_cast<Status_uint32*>    (deviceAddress + statusOffset);
	/* initialization */
	moduleAddressPtr_uint32->setValue(0);
	lengthPtr_uint32->setValue(0);
	dataInPtr_uint64[B_127_64].setValue(0);
	dataInPtr_uint64[B_63_0].setValue(0);
	dataOutPtr_uint64[B_127_64].setValue(0);
	dataOutPtr_uint64[B_63_0].setValue(0);
	/* association between the number of the register and its length */
	/*                   module,  write, read */
	ASSOCIATE_LENGTH_BIT(ADC);
	ASSOCIATE_LENGTH_BIT(DAC);
	ASSOCIATE_LENGTH_BIT(TCM);
	ASSOCIATE_LENGTH_BIT(MONITOR);
	ASSOCIATE_LENGTH_BIT(ID);
	ASSOCIATE_LENGTH_BIT(LCM1);
	ASSOCIATE_LENGTH_BIT(LCM2);
	ASSOCIATE_LENGTH_BIT(IM);
	ASSOCIATE_LENGTH_BIT(ISPM);
	ASSOCIATE_LENGTH_BIT(RESET_FAST);
}

//____________________________________________________________________________//

/**
 * Declares a register at the specified address.
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::declareRegister(P /*address*/)
{
}

/**
 * Undeclare the register at the specified address.
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::undeclareRegister(P /*address*/)
{
}

/**
 * Returns the length of the body size associated with the address for the writing
 */
template<typename T, typename P>
uint32_t AsAdControllerDevice<T, P>::lengthBodySizeOfWrite_bit()
{
	typename MemoryMap::const_iterator itReg = mem.find(reinterpret_cast<uint32_t>(moduleAddress_uint32.moduleAddress()));
	if(itReg == mem.end())
	{
		throw typename mdaq::hw::AbstractControllerDevice<T, P>::UndeclaredRegister(reinterpret_cast<P>(moduleAddress_uint32.moduleAddress()));
	}
	return itReg->second.write;
}

/**
 * Returns the length of the data associated with the address
 */
template<typename T, typename P>
uint32_t AsAdControllerDevice<T, P>::lengthBodySizeOfRead_bit()
{
	typename MemoryMap::const_iterator itReg = mem.find(reinterpret_cast<uint32_t>(moduleAddress_uint32.moduleAddress()));
	if(itReg == mem.end())
	{
		throw typename mdaq::hw::AbstractControllerDevice<T, P>::UndeclaredRegister(reinterpret_cast<P>(moduleAddress_uint32.moduleAddress()));
	}
	return itReg->second.read;
}

/**
 * Sends the write command to the state machine
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::sendCommand(void)
{
	readWritePtr_uint32->set_rw(1);
}

/**
 * Gives the local address of the register inside the module, which is the second byte
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::setLocalRegAddress()
{
	localRegAddress_uint32.setValue(inputAddress_uint32.localRegAddress());
//	cout << "localregAddress = " << localRegAddress_uint32.value() << endl;
}

/**
 * Copies the input address inside the device
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::setInputAddress(P address)
{
	inputAddress_uint32.setValue(reinterpret_cast<uint32_t>(address));
}

/**
 * Gives the address of the module, which is the first byte
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::setModuleAddress()
{
	moduleAddress_uint32.setValue(inputAddress_uint32.moduleAddress());
//	cout << "moduleAddress = " << moduleAddress_uint32.value() << endl;
}

/**
 * prints the value of the data register and of other meaningful registers
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::debug()
{
	uint32_t* dataIn32;
	uint32_t* dataOut32;

	dataIn32 = reinterpret_cast<uint32_t*>(dataInPtr_uint64);
	dataOut32 = reinterpret_cast<uint32_t*>(dataOutPtr_uint64);

	std::cout << std::noshowbase << std::hex << "======== AsAd" << std::endl;
	std::cout << "moduleAddress      = 0x" << moduleAddress_uint32.moduleAddress() << std::endl;
	std::cout << "localRegAddressPtr = 0x" << localRegAddress_uint32.value() << std::endl;
	std::cout << "lengthPtr   = " << std::dec << lengthPtr_uint32->length_bit() << std::hex << std::endl;
	std::cout << "dataIn          = 0x" << std::flush;
	for (int i = 0; i < 4; i++)
	{
		std::cout << std::hex << std::setfill('0') << std::setw(8) << *(dataIn32+i);
	}
	std::cout << std::endl << "dataOut         = 0x" << std::flush;
	for (int i = 0; i < 4; i++)
	{
		std::cout << std::hex << std::setfill('0') << std::setw(8) << *(dataOut32+i);
	}
	std::cout << std::endl << std::dec << std::setfill(' ');
}

/**
 */
template<typename T, typename P>
T AsAdControllerDevice<T, P>::readData()
{
	uint64_t data = 0;

//	cout << hex << "======== readData (AsAd)" << endl;

	switch (moduleAddress_uint32.moduleAddress())
	{
	case ADC: case LCM1 : case LCM2 : case RESET_FAST :
		throw ::mdaq::utl::CmdException("AsAd Register is Write-only");
	case DAC:
		data = dataOutPtr_uint64[B_63_0].data_DAC();
		break;
	case TCM:
		data = dataOutPtr_uint64[B_63_0].data_TCM();
		break;
	case MONITOR:
		data = dataOutPtr_uint64[B_63_0].data_MONITOR();
		break;
	case ID:
		data = dataOutPtr_uint64[B_63_0].data_ID();
		break;
	case IM:
		data = dataOutPtr_uint64[B_63_0].data_IM();
		break;
	case ISPM:
		data = dataOutPtr_uint64[B_63_0].data_ISPM();
		break;
	default:
		break;
	}
//	cout << "data = 0x" << dataOutPtr_uint64[B_63_0].secondWord() << dataOutPtr_uint64[B_63_0].firstWord() << dec << endl;
	return data;
}


/**
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>:: processBodyForReading()
{
	BodyRead_uint64 body, bodyMsb;
	body.setValue(0);
	bodyMsb.setValue(0);

	switch (moduleAddress_uint32.moduleAddress())
	{
	case ADC: case LCM1 : case LCM2 : case RESET_FAST :
		throw ::mdaq::utl::CmdException("AsAd Register is Write-only");
	case DAC:
		body.set_localRegAddress_DAC(localRegAddress_uint32.DAC());
		body.set_rw_DAC(body.reading);
		break;
	case TCM:
		body.set_localRegAddress_TCM(localRegAddress_uint32.TCM());
		body.set_rw_TCM(body.reading);
		break;
	case MONITOR:
		body.set_localRegAddress_MONITOR(localRegAddress_uint32.MONITOR());
		body.set_mode_MONITOR(body.separator);
		body.set_cmdWrite_MONITOR(body.writing_MONITOR);
		body.set_separator_MONITOR(body.separator);
		body.set_cmdRead_MONITOR(body.reading_MONITOR);
		break;
	case ID:
		body.set_localRegAddress_ID(localRegAddress_uint32.ID());
		body.set_rw_ID(body.reading);
		break;
	case IM:
		body.set_localRegAddress_IM(localRegAddress_uint32.IM());
		body.set_rw_IM(body.reading);
		break;
	case ISPM:
		body.set_localRegAddress_ISPM(localRegAddress_uint32.ISPM());
		body.set_rw_ISPM(body.reading);
		break;
	default:
		break;
	}
	dataInPtr_uint64[B_127_64].setValue(bodyMsb.value());
	dataInPtr_uint64[B_63_0].setValue(body.value());
	dataOutPtr_uint64[B_127_64].setValue(0);
	dataOutPtr_uint64[B_63_0].setValue(0);
}


/**
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::processBodyForWriting(T value)
{
	BodyWrite_uint64  body, bodyMsb;
	DataLCM_uint64 data;
	body.setValue(0);
	bodyMsb.setValue(0);

//	std::cout << "processBodyForWriting" << std::endl;
	switch (moduleAddress_uint32.moduleAddress())
	{
	case ADC:
		body.set_localRegAddress_ADC(localRegAddress_uint32.ADC());
		body.set_data_ADC(value);
		break;
	case DAC:
		body.set_localRegAddress_DAC(localRegAddress_uint32.DAC());
		body.set_data_DAC(value);
		body.set_rw_DAC(body.writing);
		break;
	case TCM:
		body.set_localRegAddress_TCM(localRegAddress_uint32.TCM());
		body.set_data_TCM(value);
		body.set_rw_TCM(body.writing);
		break;
	case MONITOR:
		body.set_localRegAddress_MONITOR(localRegAddress_uint32.MONITOR());
		body.set_data_MONITOR(value);
		body.set_mode_MONITOR(body.mode_SPI);
		body.set_cmd_MONITOR(body.writing_MONITOR);
		break;
	case ID:
		body.set_localRegAddress_ID(localRegAddress_uint32.ID());
		body.set_data_ID(value);
		body.set_rw_ID(body.writing);
		break;
	case LCM1 : case LCM2 :
		data.setValue(value);
		/* user value */
		bodyMsb.set_byte1_2_LCM(data.byte1_2());
		bodyMsb.set_msBit8_LCM(data.msBit8());
		body.set_byte3_7_LCM(data.byte3_7());
		body.set_byte10_LCM(data.byte10());
		/* constant part */
		body.set_byte8_9_LCM(body.byte8_9_LCM_c);
		body.set_conf_LCM(body.confEnable_LCM);
		break;
	case IM:
		body.set_localRegAddress_IM(localRegAddress_uint32.IM());
		body.set_data_IM(value);
		body.set_rw_IM(body.writing);
		break;
	case ISPM:
		body.set_localRegAddress_ISPM(localRegAddress_uint32.ISPM());
		body.set_data_ISPM(value);
		body.set_rw_ISPM(body.writing);
		break;
	case RESET_FAST:
		/* 0 is writen on the body */
		break;
	default:
		break;
	}
	dataInPtr_uint64[B_127_64].setValue(bodyMsb.value());
	dataInPtr_uint64[B_63_0].setValue(body.value());
	dataOutPtr_uint64[B_127_64].setValue(0);
	dataOutPtr_uint64[B_63_0].setValue(0);
}

/**
 * Waits for the busy signal of the state machine to go down
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::checkEndOfCommand()
{
	if (statusPtr_uint32->error() == 1)
	{
		throw ::mdaq::utl::CmdException("AsAd Controller Error");
	}
	else
	{
		int counter = 100;
		while ( (statusPtr_uint32->busy() == 1) && (counter > 0))
		{
			if (statusPtr_uint32->error() == 1)
			{
				throw ::mdaq::utl::CmdException("AsAd Controller Error");
			}
			counter--;
		}
		if (counter == 0)
		{
			throw ::mdaq::utl::CmdException("AsAd Controller Busy");
		}
	}
}

/**
 * Returns the value of the register at the specified address.
 */
template<typename T, typename P>
T AsAdControllerDevice<T, P>::read(P address)
{
	T value;
	setInputAddress(address);
	setModuleAddress();
	setLocalRegAddress();
	lengthPtr_uint32->set_length_bit(lengthBodySizeOfRead_bit());
	moduleAddressPtr_uint32->setValue(moduleAddress_uint32.moduleAddress());
	processBodyForReading();
	sendCommand();
	checkEndOfCommand();
	value = readData();
//	debug();

	return value;
}

/**
 * Sets the value of the register at the specified address.
 */
template<typename T, typename P>
void AsAdControllerDevice<T, P>::write(P address, T value)
{
	setInputAddress(address);
	setModuleAddress();
	setLocalRegAddress();
	lengthPtr_uint32->set_length_bit(lengthBodySizeOfWrite_bit());
	moduleAddressPtr_uint32->setValue(moduleAddress_uint32.moduleAddress());
	processBodyForWriting(value);
	sendCommand();
	checkEndOfCommand();
//	debug();
}

} //drv
} //get

#undef ASSOCIATE_LENGTH_BIT
#undef B_127_64
#undef B_63_0

#endif
