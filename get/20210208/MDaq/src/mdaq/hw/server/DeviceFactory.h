/* =============================================================================
 * DeviceFactory.h
 * -----------------------------------------------------------------------------
 * class mdaq::hw::server::DeviceFactory
 * class mdaq::hw::server::DeviceFactory::UnsupportedDevice
 * Created on: 23 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef mdaq_hw_server_DeviceFactory_h_INCLUDED
#define mdaq_hw_server_DeviceFactory_h_INCLUDED

#include "DeviceCreationParams.h"
#include "BaseDynDevice.h"
#include "../HwControllerException.h"

#include <map>
#include <string>
#include <memory>

namespace mdaq
{
namespace hw
{
namespace server
{
// forward declaration
class DeviceFactory;

typedef BaseDynDevicePtr (*DeviceFactoryMethod)(const uint64_t&);

typedef DeviceFactory* DeviceFactoryPtr;

/**
 * This class is responsible for building the adequate device objects for the
 * given register type and register width.
 */
class DeviceFactory
{
public:
	class UnsupportedDevice;

	~DeviceFactory();

	static DeviceFactory& instance();
	void registerFactoryMethod(const DeviceCreationParams& params, DeviceFactoryMethod method);
	void unregisterFactoryMethod(const DeviceCreationParams& params);
	void reset();
	BaseDynDevicePtr createDevice(const DeviceCreationParams& params, const uint64_t& deviceAddress);

private:
	DeviceFactory();

	/// a map that associate creation parameters with the factory method that
	/// instanciates the matching device.
	typedef std::map<DeviceCreationParams, DeviceFactoryMethod> FactoryMethodsMap;


	static DeviceFactoryPtr _instance;	///< the singleton unique instance.
	FactoryMethodsMap methods;
};

/**
 * This exception is thrown whenever the DeviceFactory is unable to create a
 * device matching the creation parameters.
 */
class DeviceFactory::UnsupportedDevice : public HwControllerException
{
public:
	UnsupportedDevice(const DeviceCreationParams& params);
	~UnsupportedDevice() throw();

	const DeviceCreationParams& deviceParams() const { return params; }

private:
	void writeMessage(std::ostream& stream) const;

	DeviceCreationParams params;
};

//____________________________________________________________________________//

/**
 * Returns the singleton instance, creating it if necessary.
 */
inline DeviceFactory& DeviceFactory::instance()
{
	if(_instance == 0)
	{
		_instance = new DeviceFactory;
	}
	return *_instance;
}

#define REGISTER_DEVICE_TYPE(type, policy, ctrlDev, name) \
	BaseDynDevicePtr devFactory_##type##_##ctrlDev##_##policy##_##name(const uint64_t& deviceAddress) \
	{ \
		return BaseDynDevicePtr(new DynDevice<type, policy, ctrlDev>(deviceAddress)); \
	} \
	struct DevRegistrar_##type##_##ctrlDev##_##policy##_##name \
	{ \
	    DeviceCreationParams deviceCreationParams; \
		DevRegistrar_##type##_##ctrlDev##_##policy##_##name(int w, const std::string& a, DeviceFactoryMethod m) : deviceCreationParams(w,  a)\
		{ \
			DeviceFactory::instance().registerFactoryMethod(deviceCreationParams, m); \
		} \
		~DevRegistrar_##type##_##ctrlDev##_##policy##_##name() \
		{ \
			DeviceFactory::instance().unregisterFactoryMethod(deviceCreationParams); \
		} \
	}; \
	static DevRegistrar_##type##_##ctrlDev##_##policy##_##name registrar_##type##_##ctrlDev##_##policy##_##name(sizeof(type), #name, devFactory_##type##_##ctrlDev##_##policy##_##name);

}
}
}

#endif // mdaq_hw_server_DeviceFactory_h_INCLUDED
