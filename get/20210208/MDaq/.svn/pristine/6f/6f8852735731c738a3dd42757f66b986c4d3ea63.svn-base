/* =====================================================================================================================
 * DeviceImpl.h
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
#ifndef mdaq_hw_server_DeviceImpl_h_INCLUDED
#define mdaq_hw_server_DeviceImpl_h_INCLUDED

#include "mdaq/hw/Control.h"
#include "mdaq/hw/RegisterProcessor.h"
#include "BaseDynDevice.h"

namespace mdaq
{
namespace hw
{
namespace server
{
/**
 * This class implements the ICE interface to manipulate devices.
 */
class DeviceImpl : public virtual Device
{
public:
	class Core;

private:
	BaseDynDevicePtr device;
	std::string deviceName;
	std::string regAccess;
	short regWidth;
	bool verbosity;
	void handleException();
	void completeCmdException(const std::string & regName, const std::string & fieldName="");

public:
	DeviceImpl();
	virtual ~DeviceImpl();

	std::string name(const Ice::Current & c = Ice::Current());
	Ice::Long baseAddress(const Ice::Current & c = Ice::Current());
	std::string registerAccess(const Ice::Current & c = Ice::Current());
	short registerWidth(const Ice::Current & c = Ice::Current());

	virtual void addRegister(const RegisterConfig & descr, const Ice::Current & c = Ice::Current());
	virtual void addRegisters(const RegisterConfigList & cfg, const Ice::Current & c = Ice::Current());

	Ice::Long readRegister(const std::string& registerName, const Ice::Current & c = Ice::Current());
	Ice::Long readField(const std::string& registerName, const std::string& fieldName, const Ice::Current & c = Ice::Current());

	void writeRegister(const std::string& registerName, Ice::Long value, const Ice::Current & c = Ice::Current());
	void writeField(const std::string& registerName, const std::string& fieldName, Ice::Long value, const Ice::Current & c = Ice::Current());

	void execBatch(const RegCmdSeq& cmdBatch, RegCmdSeq& batchStatus, const Ice::Current & c = Ice::Current());

	void getDeviceDescription(DeviceDescription& descr, const Ice::Current & c = Ice::Current());
	void getDeviceConfig(DeviceConfig& config, const Ice::Current & c = Ice::Current());

	void getListOfRegisters(RegisterList& registers, const Ice::Current & c = Ice::Current());
	void getListOfFields(const std::string& registerName, FieldList& fields, const Ice::Current & c = Ice::Current());

	void getRegisterConfig(const std::string& registerName, RegisterConfig& config, const Ice::Current & c = Ice::Current());
	void getRegisterDescription(const std::string& registerName, RegisterDescription & descr, const Ice::Current & c = Ice::Current());
	void getFieldDescription(const std::string& registerName, const std::string& fieldName, FieldDescription& descr, const Ice::Current & c = Ice::Current());

	void setVerbose(bool verbosity, const Ice::Current & c = Ice::Current());
	bool isVerbose(const Ice::Current & c = Ice::Current());

public:
	void init(const DeviceConfig& config);

protected:
	RegisterProcessor regProc;
};

class DeviceImpl::Core : public RegisterProcessor::Core
{
public:
	Core(DeviceImpl& deviceImpl) : deviceImpl(deviceImpl) {}
	/// @{ These are implementations of RegisterProcessor::Core interface operations (keep default implementation for other ops).
	void impl_w(const std::string& registerName, const std::string& bitFieldName, Ice::Long value);
	void impl_r(const std::string& registerName, const std::string& bitFieldName, Ice::Long& value);
	/// @}

private:
	DeviceImpl& deviceImpl;
};

}
}
}

#endif // mdaq_hw_server_DeviceImpl_h_INCLUDED
