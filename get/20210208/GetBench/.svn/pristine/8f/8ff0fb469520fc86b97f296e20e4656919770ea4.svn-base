/* =============================================================================
 * AsAdControllerDevice.h
 * -----------------------------------------------------------------------------
 * template class get::drv::AsAdControllerDevice
 * Created on: 23 September 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Joël Chavas (joel.chavas@cea.fr)
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
* complicated to manipulate, and that also ther			efore means that it is reserved
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

#ifndef get_drv_AsAdControllerDevice_h_INCLUDED
#define get_drv_AsAdControllerDevice_h_INCLUDED

#include "mdaq/hw/AbstractControllerDevice.hpp"
#include "mdaq/hw/HwControllerException.h"
#include "utl/BitFieldRegister.hpp"
#include <map>

/*================================================================
 * Defines
 *===============================================================*/

#define U_32(nb) static_cast<uint32_t>(nb)

#define DESCR(moduleName, moduleAddress, lengthWrite_bit, lengthRead_bit) \
		moduleName = U_32(moduleAddress), \
		lengthWrite_bit_##moduleName = U_32(lengthWrite_bit), \
		lengthRead_bit_##moduleName = U_32(lengthRead_bit)

#define ASSOCIATE_LENGTH_BIT(moduleName) \
{ \
	struct LengthBit_t l; \
	l.write = lengthWrite_bit_##moduleName; \
	l.read = lengthRead_bit_##moduleName; \
	mem.insert(std::make_pair(moduleName, l)); \
}

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
 *==============================================================*/

/**
 * This class implements a device policy that reads and writes from/to AsAd.
 */
template<typename T = uint64_t, typename P = T*>
class AsAdControllerDevice : public mdaq::hw::AbstractControllerDevice<T, P>
{
/*===============================================================
 * Internal offsets
 *==============================================================*/
private:
	enum Offsets {
		DEVICE_ADDRESS       = static_cast<uint32_t>(0x200000E4),
		moduleAddressOffset  = static_cast<int>(0x0),
		dataInOffset         = static_cast<int>(0x4),
		lengthOffset         = static_cast<int>(0x14),
		readWriteOffset      = static_cast<int>(0x18),
		statusOffset         = static_cast<int>(0x1C),
		dataOutOffset        = static_cast<int>(0x20)
	};

/*===============================================================
 * Constants describing the AsAd slow control protocol
 *==============================================================*/
	enum Modules {
		/*         module, address, lengthWrite_bit, lengthRead_bit */
		DESCR(        ADC,    0x00,            0x10,              0),
		DESCR(        DAC,    0x01,            0x10,           0x10),
		DESCR(        TCM,    0x02,            0x10,           0x10),
		DESCR(    MONITOR,    0x03,            0x20,           0x30),
		DESCR(         ID,    0x04,            0x10,           0x10),
		DESCR(       LCM1,    0x05,            0x51,              0),
		DESCR(       LCM2,    0x06,            0x51,              0),
		DESCR(         IM,    0x07,            0x10,           0x10),
		DESCR(       ISPM,    0x08,            0x10,           0x10),
		DESCR( RESET_FAST,    0x0F,            0x10,           0x10)
	}; ///< list of the AsAd modules identified by their address

/*===============================================================
 * Definitions of the bitfields
 *==============================================================*/

	typedef utl::BitFieldRegister<volatile uint32_t> BitField_uint32;
	typedef utl::BitFieldRegister<volatile uint64_t> BitField_uint64;

	/**
	 */
	struct Status_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(checksumError,    3, 1)
		DECLARE_FIELD(error,            2, 1)
		DECLARE_FIELD(busy,             1, 1)
		DECLARE_FIELD(valid,            0, 1)
	};

	struct BitFieldDescr_uint64 : public BitField_uint64
	{
		DECLARE_FIELD(secondWord,  32, 32)
		DECLARE_FIELD(firstWord,    0, 32)
	};

	/**
	 * Defines the length of the register address for each module
	 */
	struct LocalRegAddress_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(ADC,           0, 5)
		DECLARE_FIELD(DAC,           0, 5)
		DECLARE_FIELD(TCM,           0, 7)
		DECLARE_FIELD(MONITOR,       0, 8)
		DECLARE_FIELD(ID,            0, 7)
		DECLARE_FIELD(IM,            0, 7)
		DECLARE_FIELD(ISPM,          0, 7)
		DECLARE_FIELD(RESET_FAST,    0, 7)
	};

	/**
	 * Defines the position inside the body of the slow control protocol
	 * of the local register address, of the data and of the read/write bit
	 * for the writing
	 */
	struct BodyWrite_uint64 : public BitField_uint64
	{
		static const uint64_t reading = 1;
		static const uint64_t writing = 0;
		static const uint64_t mode_I2C = 0X00;
		static const uint64_t mode_SPI = 0X80;
		static const uint64_t writing_MONITOR = 0x90;
		static const uint64_t byte8_9_LCM_c = 0x0001;
		static const uint64_t confEnable_LCM = 1;
	/* localRegAddress */
		DECLARE_FIELD(localRegAddress_ADC,  11, 5)
		DECLARE_FIELD(localRegAddress_DAC,  8, 5)
		DECLARE_FIELD(localRegAddress_TCM,  8, 7)
		DECLARE_FIELD(localRegAddress_MONITOR,  8, 8) /* write */
		DECLARE_FIELD(localRegAddress_ID,   8, 7)
		DECLARE_FIELD(localRegAddress_IM,   8, 7)
		DECLARE_FIELD(localRegAddress_ISPM, 8, 7)
	/* data */
		DECLARE_FIELD(data_ADC,  0, 11)
		DECLARE_FIELD(data_DAC,  0, 8)
		DECLARE_FIELD(data_TCM,  0, 8)
		DECLARE_FIELD(data_MONITOR,  0, 8) /* write */
		DECLARE_FIELD(data_ID,   0, 8)
		DECLARE_FIELD(data_IM,   0, 8)
		DECLARE_FIELD(data_ISPM, 0, 8)
	/* configuration MONITOR */
		DECLARE_FIELD(mode_MONITOR, 24, 8);
		DECLARE_FIELD(cmd_MONITOR, 16, 8);
	/* configuration LCM */
		DECLARE_FIELD(byte1_2_LCM, 1, 16); /* LSB word, at position 64+1 in the body */
		DECLARE_FIELD(msBit8_LCM, 0, 1); /* LSB word, at position 64+0 in the body */
		DECLARE_FIELD(byte3_7_LCM, 25, 39);
		DECLARE_FIELD(byte8_9_LCM, 9, 16);
		DECLARE_FIELD(byte10_LCM, 1, 8);
		DECLARE_FIELD(conf_LCM, 0, 1);
	/* rw bits */
		DECLARE_FIELD(rw_DAC,  15, 1)
		DECLARE_FIELD(rw_TCM,  15, 1)
		DECLARE_FIELD(rw_ID,   15, 1)
		DECLARE_FIELD(rw_IM,   15, 1)
		DECLARE_FIELD(rw_ISPM, 15, 1)
	};

	/**
	 * Defines the position inside the body of the slow control protocol
	 * of the local register address, of the data and of the read/write bit
	 * for the reading
	 */
	struct BodyRead_uint64 : public BitField_uint64
	{
		static const uint64_t reading = 1;
		static const uint64_t writing = 0;
		static const uint64_t reading_MONITOR = 0x91;
		static const uint64_t writing_MONITOR = 0x90;
		static const uint64_t separator = 0xFF;
		DECLARE_FIELD(secondWord,  32, 32)
		DECLARE_FIELD(firstWord,  0, 32)
	/* localRegAddress */
		DECLARE_FIELD(localRegAddress_DAC,  8, 5)
		DECLARE_FIELD(localRegAddress_TCM,  8, 7)
		DECLARE_FIELD(localRegAddress_MONITOR,  24, 8) /* read */
		DECLARE_FIELD(localRegAddress_ID,   8, 7)
		DECLARE_FIELD(localRegAddress_IM,   8, 7)
		DECLARE_FIELD(localRegAddress_ISPM, 8, 7)
	/* data */
		DECLARE_FIELD(data_DAC,  0, 8)
		DECLARE_FIELD(data_TCM,  0, 8)
		DECLARE_FIELD(data_MONITOR,  0, 8) /* read */
		DECLARE_FIELD(data_ID,   0, 8)
		DECLARE_FIELD(data_IM,   0, 8)
		DECLARE_FIELD(data_ISPM, 0, 8)
	/* configuration MONITOR */
		DECLARE_FIELD(mode_MONITOR, 40, 8);
		DECLARE_FIELD(cmdWrite_MONITOR, 32, 8);
		DECLARE_FIELD(separator_MONITOR, 16, 8);
		DECLARE_FIELD(cmdRead_MONITOR, 8, 8);
	/* rw bits */
		DECLARE_FIELD(rw_DAC,  15, 1)
		DECLARE_FIELD(rw_TCM,  15, 1)
		DECLARE_FIELD(rw_ID,   15, 1)
		DECLARE_FIELD(rw_IM,   15, 1)
		DECLARE_FIELD(rw_ISPM, 15, 1)
	};

	/**
	 * Splits the 64-bit input data
	 * among the two 64_bit field of the AsAd slow control body
	 */
	struct DataLCM_uint64 : public BitField_uint64
	{
		DECLARE_FIELD(byte1_2, 48, 16)
		DECLARE_FIELD(msBit8 , 47, 1)
		DECLARE_FIELD(byte3_7,  8, 39)
		DECLARE_FIELD(byte10,   0, 8)
	};

	struct InputAddress_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(moduleAddress,    8, 8)
		DECLARE_FIELD(localRegAddress,  0, 8)
	};

	struct ModuleAddress_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(moduleAddress,  0, 8)
	};

	struct ReadWrite_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(rw,  0, 1)
	};

	/**
	 * Describes the length of the body size in bits
	 * It is at most 0x51 wide. 7 bits are enough to memorize it
	 */
	struct Length_uint32 : public BitField_uint32
	{
		DECLARE_FIELD(length_bit,  0, 7)
	};

/*================================================================
 * Misc
 *===============================================================*/

	struct LengthBit_t
	{
		uint32_t write;
		uint32_t read;
	};

/*===============================================================
 * Device address
 *==============================================================*/
	char* deviceAddress;

/*===============================================================
 * Public functions
 *==============================================================*/
public:

	typedef std::map<uint32_t, struct LengthBit_t> MemoryMap;

	MemoryMap mem;

	AsAdControllerDevice(uint64_t &baseAddress);
	virtual ~AsAdControllerDevice() {}

	void declareRegister(P address);
	void undeclareRegister(P address);
	uint32_t lengthBodySizeOfRead_bit();
	uint32_t lengthBodySizeOfWrite_bit();
	void sendCommand(void);
	void checkEndOfCommand(void);
	void setInputAddress(P address);
	void setModuleAddress();
	void setLocalRegAddress();
	void processBodyForReading();
	void processBodyForWriting(T value);
	T readData();
	T read(P address);
	void write(P address, T value);
	void debug();

/*===============================================================
 * Pointers to firmware registers
 *==============================================================*/
private:
	Length_uint32*          lengthPtr_uint32;
	ReadWrite_uint32*       readWritePtr_uint32;
	Status_uint32*          statusPtr_uint32;
	ModuleAddress_uint32*   moduleAddressPtr_uint32;
	BodyWrite_uint64*       dataInPtr_uint64;
	BodyRead_uint64*        dataOutPtr_uint64;
	/* local variables */
	InputAddress_uint32     inputAddress_uint32;
	LocalRegAddress_uint32  localRegAddress_uint32;
	ModuleAddress_uint32    moduleAddress_uint32;
};

}
}

#undef U_32
#undef DESCR

#endif
