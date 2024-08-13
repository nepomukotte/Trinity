/* =====================================================================================================================
 * HwNodeImpl.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::hw::server::HwNodeImpl
 * Created on: 19 August 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of
 * the Mordicus Real-Time Software Development Framework.
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

#include "mdaq/hw/server/HwNodeImpl.h"
#include "mdaq/hw/server/DeviceImpl.h"
#include "utl/Logging.h"
#include ICE_HEADER(ObjectAdapter.h)
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace mdaq
{
namespace hw
{
namespace server
{
/**
 * Constructor.
 */
HwNodeImpl::HwNodeImpl() : nodeName()
{
}

/**
 * Destructor.
 */
HwNodeImpl::~HwNodeImpl()
{
}

/**
 * Says hello for debug
 */
void HwNodeImpl::testConnectionToHardware(const Ice::Current&)
{
	std::cout << "Hello!" << std::endl;
}

/**
 * Shuts down the communicator, and thus closes the server
 */
void HwNodeImpl::shutdown(const Ice::Current& c)
{
	LOG_INFO() << "Shutting down...";
    c.adapter->getCommunicator()->shutdown();
}

/**
 * Reboots the board  (VXWORKS) or shuts down the communicator (other systems)
 */
void HwNodeImpl::reboot(const Ice::Current & c)
{
	LOG_INFO() << "Rebooting...";
#ifdef VXWORKS
    IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(1000));
    uint32_t rebootFuncPtr=0xFFFFFFFC;
    ((void (*)(void))rebootFuncPtr)();
#else
    shutdown(c);
#endif
}

/**
 * Returns the name of the node.
 */
std::string HwNodeImpl::name(const Ice::Current& current)
{
	return nodeName.empty() ? current.adapter->getName() : nodeName;
}

/**
 * Sets the name of the node.
 * @param name Name to give to the node.
 */
void HwNodeImpl::setName(const ::std::string & name, const ::Ice::Current&)
{
	static int nameCount = 0;
	nodeName = name;
	LOG_INFO() << std::string(51, '=') << ' ' << name << '-' << nameCount++ << ' ' << std::string(51, '=');
}

/**
 * creates the node initializing with parameters describing devices, registers and fields.
 */
void HwNodeImpl::create(const NodeConfig& config, const Ice::Current& current)
{
	try
	{
		nodeName = config.id;
		LOG_DEBUG() << "Creating node: " << nodeName;
		adapter_ = current.adapter;
		for(DeviceConfigList::const_iterator it = config.devices.begin(),
			end = config.devices.end(); it != end; ++it)
		{
			createDevice(*it, current);
		}
	}
	catch(HwControllerException& e)
	{
		throw ::mdaq::utl::CmdException(e.getMessage());
	}
	catch (const std::exception & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
	catch (...)
	{
		throw ::mdaq::utl::CmdException("Unknown exception while creating node!");
	}
}

/**
 * Adds a device and its registers to the node with given device configuration.
 * @param config Configuration of the device.
 */
void HwNodeImpl::createDevice(const DeviceConfig & config, const Ice::Current& current)
{
	// LOG_DEBUG() << "Creating device: " << config.descr.name;
	DeviceImpl* devImpl = new DeviceImpl;
	devImpl->init(config);
	DevicePrx deviceProxy = DevicePrx::uncheckedCast(current.adapter->add(devImpl, makeDeviceIdentity(config.descr.name)));
	deviceMap[config.descr.name] = deviceProxy;
	// LOG_DEBUG() << "Device " << config.descr.name << " created";
}

/**
 * Adds a device to the node with given device configuration.
 * @param deviceDescr Description of the device (name, policy, base address and register width)
 */
void HwNodeImpl::addDevice(const mdaq::hw::DeviceDescription & deviceDescr, const Ice::Current & c)
{
	try
	{
		// Remove device with same name if it exists
		removeDevice(deviceDescr.name);

		// Create device configuration (description plus empty list of registers)
		DeviceConfig cfg;
		cfg.descr = deviceDescr;

		// Create device
		createDevice(cfg, c);
	}
	catch(HwControllerException& e)
	{
		throw ::mdaq::utl::CmdException(e.getMessage());
	}
	catch (const std::exception & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
}

/**
 * Destroys the node.
 */
void HwNodeImpl::destroy(const Ice::Current& current)
{
	removeAllDevices(current);
}

/**
 * Return a map of the proxies to devices within this node, indexed by their name.
 */
void HwNodeImpl::getMapOfDevices(DeviceMap & map, const ::Ice::Current&)
{
	map = deviceMap;
}

/**
 * Returns the names of devices within this node.
 */
void HwNodeImpl::getListOfDevices(DeviceList & list, const Ice::Current&)
{
	list.clear();
	for (DeviceMap::iterator it = deviceMap.begin(); it != deviceMap.end(); ++it)
	{
		list.push_back(it->first);
	}
}

void HwNodeImpl::getNodeDescription(std::string & nodeId, const Ice::Current& current)
{
	nodeId = name(current);
}

void HwNodeImpl::getNodeConfig(NodeConfig &, const Ice::Current&)
{
	throw ::mdaq::utl::CmdException("This feature has yet to be implemented.");
}

/**
 * Returns the proxy of the device having the specified name.
 */
DevicePrx HwNodeImpl::findDevice(const std::string & name, const Ice::Current &)
{
	// LOG_DEBUG() << "Looking for device: " << name;
	DeviceMap::const_iterator deviceIt = deviceMap.find(name);
	if(deviceIt == deviceMap.end())
	{
		LOG_ERROR() << "Device not found: " << name;
		throw ::mdaq::utl::CmdException(std::string("Device not found: ") + name);
	}
	// LOG_DEBUG() << "Returning proxy to device '" << name << "' of node '" << nodeName << '\'';
	return deviceIt->second;
}

/**
 * Returns the local device having the specified name.
 */
DeviceImpl & HwNodeImpl::findLocalDevice(const std::string & name)
{
	// LOG_DEBUG() << "Looking for local device: " << name;
	if (not adapter_)
	{
		LOG_ERROR() << "Device not found: " << name;
		throw ::mdaq::utl::CmdException(std::string("Device not found: ") + name);
	}
	Ice::ObjectPtr deviceObj = adapter_->find(makeDeviceIdentity(name));
	if (not deviceObj)
	{
		LOG_ERROR() << "Device not found: " << name;
		throw ::mdaq::utl::CmdException(std::string("Device not found: ") + name);
	}
	return dynamic_cast< DeviceImpl& >(*deviceObj);
}

/**
 * Removes the device with given name from the node.
 * @param deviceName Name of the device to remove.
 */
void HwNodeImpl::removeDevice(const std::string & deviceName, const Ice::Current & current)
{
	DeviceMap::iterator device = deviceMap.find(deviceName);
	if (device != deviceMap.end())
	{
		try
		{
			current.adapter->remove(makeDeviceIdentity(device->first));
		}
		catch (const Ice::NotRegisteredException & e)
		{
			throw ::mdaq::utl::CmdException(
			        std::string("Error removing device '") + device->first + "': " + e.ice_name());
		}
		deviceMap.erase(device);
	}
}

/**
 * Removes all devices from the node.
 */
void HwNodeImpl::removeAllDevices(const Ice::Current& current)
{
	for (DeviceMap::const_iterator it = deviceMap.begin(), end = deviceMap.end();
		it != end; ++it)
	{
		try
		{
			current.adapter->remove(makeDeviceIdentity(it->first));
		}
		catch (const Ice::NotRegisteredException & e)
		{
			throw ::mdaq::utl::CmdException(std::string("Error removing device '") + it->first + "': " + e.ice_name());
		}
	}
	deviceMap.clear();
}

/**
 * Creates an ICE identity for a device having the specified name.
 */
Ice::Identity HwNodeImpl::makeDeviceIdentity(const std::string& deviceName)
{
	Ice::Identity id;
	id.name  = deviceName;
	return id;
}

} // namespace server
} // namespace hw
} // namespace mdaq
