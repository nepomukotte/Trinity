/* =====================================================================================================================
 * DeviceImpl.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::hw::server::DeviceImpl
 * Created on: 19 August 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
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

#include "DeviceImpl.h"
#include "DeviceFactory.h"
#include "mdaq/utl/CmdException.h"
#include "utl/Logging.h"

#include ICEUTIL_HEADER(Thread.h)
#include ICEUTIL_HEADER(Time.h)
#include ICE_HEADER(LocalException.h)
#include <sstream>

// Macro to handle exception (exception dispatcher pattern)
#define HANDLE_LOOKUP_EXCEP(instruction)		\
	try											\
	{											\
		instruction;							\
	}											\
	catch (...) 								\
	{											\
		handleException();						\
	}

#define COMPLETE_CMD_EXCEP(cmd, regName, fieldName)	\
	try												\
	{												\
		cmd;										\
	}												\
	catch (...) 									\
	{												\
		completeCmdException(regName, fieldName);	\
	}

#define DBG_IF(cond) if(UTL_LOG_LEVEL > 0 or !(cond)) ; else ::utl::BackendLogger(::utl::debug, __FILE__, __LINE__)

using ::mdaq::hw::RegisterProcessor;

namespace mdaq
{
namespace hw
{
namespace server
{
/**
 * Constructor.
 */
DeviceImpl::DeviceImpl() : regWidth(0), verbosity(false), regProc(RegisterProcessor::CorePtr(new DeviceImpl::Core(*this)))
{
}

/**
 * Destructor.
 */
DeviceImpl::~DeviceImpl()
{
}

/**
 * Re-throws exception after adding device, register and field names to error message.
 */
/*void DeviceImpl::handleRegisterCmdException(const std::string & registerName, const std::string fieldName)
{
	try
	{
		throw;
	}
	catch (const ::mdaq::utl::CmdException & e)
	{
		std::ostringstream oss;
		oss << "Error accessing ";
		if (not fieldName.empty())
		{
			oss << "field '" << fieldName << "' of ";
		}
		oss << "register '" << registerName << "' of device '" << deviceName << "'(" << e.reason << ')';
		throw ::mdaq::utl::CmdException(oss.str());
	}
}*/

/**
 * Returns the name of the device.
 */
std::string DeviceImpl::name(const Ice::Current&)
{
	return deviceName;
}

/**
 * Returns the base address of the device.
 */
Ice::Long DeviceImpl::baseAddress(const Ice::Current&)
{
	return device->baseAddress();
}

/**
 * Returns the title associated with the type of access of the registers.
 */
std::string DeviceImpl::registerAccess(const Ice::Current&)
{
	return regAccess;
}

/**
 * Returns the size of the registers of the device.
 */
short DeviceImpl::registerWidth(const Ice::Current&)
{
	return regWidth;
}

/**
 * Returns the value of the specified register.
 */
Ice::Long DeviceImpl::readRegister(const std::string& registerName, const Ice::Current&)
{
	Ice::Long value = -1;
	HANDLE_LOOKUP_EXCEP(COMPLETE_CMD_EXCEP(regProc.read(registerName, "", value), registerName, ""))
	return value;
}

/**
 * Returns the value of the specified field.
 */
Ice::Long DeviceImpl::readField(const std::string& registerName, const std::string& fieldName, const Ice::Current&)
{
	Ice::Long value = -1;
	HANDLE_LOOKUP_EXCEP(COMPLETE_CMD_EXCEP(regProc.read(registerName, fieldName, value), registerName, fieldName))
	return value;
}

/**
 * Changes the value of the specified register.
 */
void DeviceImpl::writeRegister(const std::string& registerName, Ice::Long value, const Ice::Current&)
{
	HANDLE_LOOKUP_EXCEP(COMPLETE_CMD_EXCEP(regProc.write(registerName, "", value), registerName, ""))
}

/**
 * Changes the value of the specified field.
 */
void DeviceImpl::writeField(const std::string& registerName, const std::string& fieldName, Ice::Long value, const Ice::Current&)
{
	HANDLE_LOOKUP_EXCEP(COMPLETE_CMD_EXCEP(regProc.write(registerName, fieldName, value), registerName, fieldName))
}

/**
 * Executes a batch of register commands for improved execution performance.
 * @param cmdBatch Vector of read/write commands to be executed.
 * @param batchStatus Vector of execution status of batch commands (after their execution).
 * @return Execution results: for read commands the result of the read in "parameter" field, in both read and write
 *         case, the exception status in "status" field.
 */
void DeviceImpl::execBatch(const RegCmdSeq& cmdBatch, RegCmdSeq& batchStatus, const Ice::Current&)
{
	HANDLE_LOOKUP_EXCEP(regProc.exec(cmdBatch, batchStatus))
}

/**
 * Enumerates the device's register by making a list with the register names.
 */
void DeviceImpl::getListOfRegisters(RegisterList& registers, const Ice::Current&)
{
	HANDLE_LOOKUP_EXCEP(device->registerList(registers))
}

/**
 * Enumerates a registers' fields by making a list with the field names.
 */
void DeviceImpl::getListOfFields(const std::string& registerName, FieldList& fields, const Ice::Current&)
{
	HANDLE_LOOKUP_EXCEP(device->registerFind(registerName).fieldList(fields))
}

/**
 * Returns the description of the specified register.
 * @param registerName Name of the register.
 * @param,out descr Description of the register with specified name.
 */
void DeviceImpl::getRegisterDescription(const std::string& registerName, RegisterDescription & descr, const Ice::Current&)
{
	descr.name = registerName;
	HANDLE_LOOKUP_EXCEP(descr.offset = device->registerFind(registerName).storageAddress() - device->baseAddress())
}

/**
 * Returns the description of the specified field in the specified register.
 * @param registerName Name of the register.
 * @param fieldName Name of the field.
 * @param,out descr Description of the field.
 */
void DeviceImpl::getFieldDescription(const std::string& registerName, const std::string& fieldName, FieldDescription& descr, const Ice::Current&)
{
	descr.name = fieldName;
	HANDLE_LOOKUP_EXCEP(descr.offset = device->registerFind(registerName).fieldOffset(fieldName))
	HANDLE_LOOKUP_EXCEP(descr.width = device->registerFind(registerName).fieldWidth(fieldName))
}

/**
 * Initializes registers and fields from the specified configuration data.
 */
void DeviceImpl::init(const DeviceConfig& config)
{
	deviceName = config.descr.name;
	regAccess = config.descr.registerAccess;
	regWidth = config.descr.registerWidth;

	{	// hack for vxworks which doesn't support assignment from return!!
		BaseDynDevicePtr tmpDevice = DeviceFactory::instance().createDevice(
			DeviceCreationParams(config.descr.registerWidth, config.descr.registerAccess), config.descr.baseAddress);
		device = tmpDevice;
	}

	// Add registers and their fields
	addRegisters(config.registers);
}

/**
 * Adds a register and its fields
 * @param cfg register description (name, offset, fields)
 */
void DeviceImpl::addRegister(const RegisterConfig & cfg, const Ice::Current &)
{
	HANDLE_LOOKUP_EXCEP(
		BaseDynRegister& reg = device->registerCreate(cfg.descr.name, cfg.descr.offset);

		for (FieldDescriptionList::const_iterator itField = cfg.fields.begin(),
			endField = cfg.fields.end(); itField != endField; ++itField)
		{
			reg.fieldAdd(itField->name, itField->offset, itField->width);
		}
	)
}

/**
 * Adds a list of registers and its fields
 * @param cfg list of register descriptions (name, offset, fields)
 */
void DeviceImpl::addRegisters(const RegisterConfigList & cfg, const Ice::Current &)
{
	// Add registers and their fields
	for (RegisterConfigList::const_iterator itReg = cfg.begin(),
		endReg = cfg.end(); itReg != endReg; ++itReg)
	{
		addRegister(*itReg);
	}
}

/**
 * Returns the description of the device.
 * @param[out] descr Device description.
 */
void DeviceImpl::getDeviceDescription(DeviceDescription & descr, const Ice::Current &)
{
	HANDLE_LOOKUP_EXCEP (
		descr.name = this->deviceName;
		descr.baseAddress = device->baseAddress();
		descr.registerAccess = this->regAccess;
		descr.registerWidth = this->regWidth;
	)
}

/**
 * Returns the configuration of the device.
 * @param[out] config Device configuration.
 */
void DeviceImpl::getDeviceConfig(DeviceConfig & config, const Ice::Current & current)
{
	getDeviceDescription(config.descr, current);

	RegisterList registerList;
	getListOfRegisters(registerList, current);

	config.registers.resize(registerList.size());
	for (size_t i=0; i < registerList.size(); ++i)
	{
		getRegisterConfig(registerList[i], config.registers[i], current);
	}

}

void DeviceImpl::getRegisterConfig(const std::string & registerName, RegisterConfig & config, const Ice::Current & current)
{
	getRegisterDescription(registerName, config.descr, current);

	FieldList fieldList;
	getListOfFields(registerName, fieldList, current);

	config.fields.resize(fieldList.size());
	for (size_t i=0; i < fieldList.size(); ++i)
	{
		getFieldDescription(registerName, fieldList[i], config.fields[i], current);
	}
}

void DeviceImpl::setVerbose(bool verbosity, const Ice::Current &)
{
	this->verbosity = verbosity;
}

bool DeviceImpl::isVerbose(const Ice::Current&)
{
	return verbosity;
}

void DeviceImpl::handleException()
{
	try
	{
		throw;
	}
	catch (const BaseDynDevice::RegisterNotFound& e)
	{
		LOG_ERROR() << "RegisterNotFound: " << e.getMessage();
		throw ::mdaq::utl::CmdException(e.getMessage());
	}
	catch (const BaseDynRegister::FieldNotFound& e)
	{
		LOG_ERROR() << "FieldNotFound: " << e.getMessage();
		throw ::mdaq::utl::CmdException(e.getMessage());
	}
	catch (const BaseDynRegister::OutOfRange& e)
	{
		LOG_ERROR() << "OutOfRange: " << e.getMessage();
		throw ::mdaq::utl::CmdException(e.getMessage());
	}
	catch (const Ice::ObjectNotExistException & e)
	{
		LOG_ERROR() << "Device '" << e.id.name << "' does not exist";
		throw ::mdaq::utl::CmdException(std::string("Device '") + e.id.name + "' does not exist");
	}
	catch (const ::mdaq::utl::CmdException& e)
	{
		LOG_ERROR() << "CmdException: " << e.reason;
		throw;
	}
	catch (const std::exception& e)
	{
		LOG_ERROR() << "Exception: " << e.what();
		throw ::mdaq::utl::CmdException(e.what());
	}
	catch (const std::string& msg)
	{
		LOG_ERROR() << msg;
		throw ::mdaq::utl::CmdException(msg);
	}
	catch (const char* msg)
	{
		LOG_ERROR() << msg;
		throw ::mdaq::utl::CmdException(msg);
	}
	catch (...)
	{
		LOG_ERROR() << "Unknown Exception";
		throw ::mdaq::utl::CmdException("Unknown Exception");
	}
}

/**
 * Completes command error with name of device, register and field.
 * Uses the exception dispatcher pattern.
 */
void DeviceImpl::completeCmdException(const std::string & regName, const std::string & fieldName)
{
	try
	{
		throw;
	}
	catch (const ::mdaq::utl::CmdException & e)
	{
		std::ostringstream oss;
		oss << "Error accessing ";
		if (not fieldName.empty())
		{
			oss << "field '" << fieldName << "' of ";
		}
		oss << "register '" << regName << "' of device '" << deviceName << "' (" << e.reason << ')';
		throw ::mdaq::utl::CmdException(oss.str());
	}
}

void DeviceImpl::Core::impl_w(const std::string& registerName, const std::string& fieldName, Ice::Long value)
{
	DBG_IF(deviceImpl.verbosity) << "writeField(" << registerName << "," << fieldName  << ", "<< std::hex << value << ")" << std::dec;
	if (not fieldName.empty())
	{
		deviceImpl.device->registerFind(registerName).fieldWriteInt64(fieldName, value);
	}
	else
	{
		deviceImpl.device->registerFind(registerName).setValueInt64(value);
	}
}

void DeviceImpl::Core::impl_r(const std::string& registerName, const std::string& fieldName, Ice::Long& value)
{
	DBG_IF(deviceImpl.verbosity) << "readField(" << registerName << "," << fieldName << ")";
	if (not fieldName.empty())
	{
		value = deviceImpl.device->registerFind(registerName).fieldReadInt64(fieldName);
	}
	else
	{
		value = deviceImpl.device->registerFind(registerName).getValueInt64();
	}
}

}
}
}
