/* =====================================================================================================================
* AGetControllerDevice.h
 * ---------------------------------------------------------------------------------------------------------------------
* spacialized template class get::drv::AGetControllerDevice<AGetWord, AGetAddr>
* Created on: 23 September 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of the Mordicus Real-Time Software
 * Development Framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef get_drv_AGetControllerDevice_h_INCLUDED
#define get_drv_AGetControllerDevice_h_INCLUDED

#include "mdaq/hw/AbstractControllerDevice.hpp"
#include "mdaq/hw/HwControllerException.h"
#include "utl/BitFieldRegister.hpp"

/*================================================================
 * Defines
 *===============================================================*/

#define U_32(nb) static_cast<uint32_t>(nb)

#define DESCR(regName, regAddress, length_bit) \
		regName = U_32(regAddress), \
		length_##regName = U_32(length_bit), \
		offset_##regName = U_32(64-length_bit)

#define ASSOCIATE_LENGTH_BIT(regName) \
	mem.insert(std::make_pair(regName, length_##regName))

#define B_127_64 0
#define B_63_0 1

/*================================================================
 * Namespace get::drv
 *===============================================================*/

namespace get
{
namespace drv
{

/*===============================================================
* Class declaration
*===============================================================*/

/**
 * This class implements a device policy that simulates reading and writing
 * registers into an addressed space.
 */
template<typename T = uint64_t, typename P = T*>
class AGetControllerDevice : public mdaq::hw::AbstractControllerDevice<T, P>
{
/*===============================================================
 * Internal offsets
 *==============================================================*/
private:
	enum Offsets {
		// Offsets from address 0x20000000
		regAddressOffset = static_cast<int>(0xA0),  // Nathan: agetRegisterAddress
		dataInOffset     = static_cast<int>(0xA4),  // Nathan: dataToAget3 (corresponding to AGET0 on AsAd?)
		lengthOffset     = static_cast<int>(0xB4), // Nathan: aget_sc_size
		readWriteOffset  = static_cast<int>(0xB8), // Nathan: aget_sc_rw
		dataOutOffset    = static_cast<int>(0xC0), // Nathan: dataFromAget3 (corresponding to AGET0 on AsAd?)
		statusOffset     = static_cast<int>(0xBC)  // Nathan: agetControllerStatus

		// Nathan: agetChipSelect        @ address 0x20000070
		// Nathan: aget_sc_configuration @ address 0x20000074
	};


/*===============================================================
 * Constants describing the AGet slow control protocol
 *==============================================================*/
	enum Modules {
		/*    register, address, length_bit */
		DESCR(    reg0,       0,         64),
		DESCR(    reg1,       1,         32),
		DESCR(    reg2,       2,         32),
		DESCR(    reg3,       3,         34),
		DESCR(    reg4,       4,         34),
		DESCR(    reg5,       5,         16),
		DESCR(    reg6,       6,         64),
		DESCR(    reg7,       7,         64),
		DESCR(    reg8,       8,        128),
		DESCR( reg8msb,      13,        128),
		DESCR(    reg9,       9,        128),
		DESCR( reg9msb,      14,        128),
		DESCR(   reg10,      10,         64),
		DESCR(   reg11,      11,         64),
		DESCR(   reg12,      12,         16),
		DESCR(     hit,      15,         68),
		DESCR(  hitmsb,      16,         68),
	}; ///< list of the AsAd modules identified by their address

/*===============================================================
 * Definitions of the bitfields
 *==============================================================*/

	typedef utl::BitFieldRegister<volatile uint32_t> BitField_uint32;
	typedef utl::BitFieldRegister<volatile uint64_t> BitField_uint64;

	/**
	 * BitField structure describing AGetControllerStatus
	 */
	struct Status_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(error            , 2, 1)
		DECLARE_FIELD(busy             , 1, 1)
		DECLARE_FIELD(valid            , 0, 1)
	};

	/**
	 * Defines the position inside the body of the slow control protocol
	 * of the local register address, of the data and of the read/write bit
	 * for the writing
	 */
	struct BodyWrite_uint64 : public BitField_uint64
	{
		DECLARE_FIELD( reg0   ,     offset_reg0,  length_reg0)
		DECLARE_FIELD( reg1   ,     offset_reg1,  length_reg1)
		DECLARE_FIELD( reg2   ,     offset_reg2,  length_reg2)
		DECLARE_FIELD( reg3   ,     offset_reg3,  length_reg3)
		DECLARE_FIELD( reg4   ,     offset_reg4,  length_reg4)
		DECLARE_FIELD( reg5   ,     offset_reg5,  length_reg5)
		DECLARE_FIELD( reg6   ,     offset_reg6,  length_reg6)
		DECLARE_FIELD( reg7   ,     offset_reg7,  length_reg7)
		DECLARE_FIELD( reg8   ,               0,           64)
		DECLARE_FIELD( reg8msb,               0,           64)
		DECLARE_FIELD( reg9   ,               0,           64)
		DECLARE_FIELD( reg9msb,               0,           64)
		DECLARE_FIELD( reg10  ,    offset_reg10, length_reg10)
		DECLARE_FIELD( reg11  ,    offset_reg11, length_reg11)
		DECLARE_FIELD( reg12  ,    offset_reg12, length_reg12)
		DECLARE_FIELD( msb4   ,              60,            4)
		DECLARE_FIELD( lsb60  ,               0,           60)
	};

	/**
	 * Defines the position inside the body of the slow control protocol
	 * of the local register address, of the data and of the read/write bit
	 * for the reading
	 */
	struct BodyRead_uint64 : public BitField_uint64
	{
		DECLARE_FIELD( lsb4   ,              0,            4)
		DECLARE_FIELD( msb60  ,              4,           60)
	};

	struct RegAddress_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(regAddress,  0, 7)
	};

	struct Length_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(regAddress,  0, 7)
	};

	struct ReadWrite_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(rw,  0, 1)
	};



/*===============================================================
* Device address
*==============================================================*/
	char* deviceAddress;

/*===============================================================
 * Public functions
 *==============================================================*/
public:
	typedef std::map<uint32_t, uint32_t> MemoryMap;

	MemoryMap mem;

	AGetControllerDevice(uint64_t & baseAddress);
	virtual ~AGetControllerDevice() {}

	void declareRegister(P address);
	void undeclareRegister(P address);
	uint32_t length(P address);
	void checkIfReadingIsAllowed(P address);
	void checkIfWritingIsAllowed(P address);
	void setLength(P address);
	void setRegAddress(P address);
	void assertNotBusy(size_t numTrials = 10000);
	void checkEndOfCommand();
	void processBodyForReading();
	T processBodyAfterReading();
	void processBodyForWriting(T value);
	void sendReadCommand();
	void sendWriteCommand();
	void debug();;
	T read(P address);
	void write(P address, T value);

/*===============================================================
 * Pointers to firmware registers
 *==============================================================*/
private:
	RegAddress_uint32* regAddressPtr_uint32;
	Length_uint32*     lengthPtr_uint32;
	ReadWrite_uint32*  readWritePtr_uint32;
	BodyWrite_uint64*  dataInPtr_uint64;
	BodyRead_uint64*   dataOutPtr_uint64;
	Status_uint32*     statusPtr_uint32;
	/* local variable */
	uint32_t           fictiveAddress_uint32;
};

} //drv
} //get

#undef DESCR
#undef U_32


#endif //get_drv_AGetControllerDevice_h_INCLUDED
