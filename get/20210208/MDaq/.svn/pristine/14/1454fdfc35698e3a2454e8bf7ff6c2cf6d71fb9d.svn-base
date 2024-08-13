/* =============================================================================
 * DeviceFactory.cpp
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

#include "DeviceFactory.h"
#include <utility>

namespace mdaq
{
namespace hw
{
namespace server
{
/**
 * The singleton instance.
 */
DeviceFactoryPtr DeviceFactory::_instance = 0;

/**
 * Constructor.
 */
DeviceFactory::DeviceFactory()
{
}

/**
 * Destructor.
 */
DeviceFactory::~DeviceFactory()
{
}

/**
 * Registers a method capable of creating a device supporting the specified
 * creation parameters.
 * If the factory method already exist for these parameters, it will be
 * replaced by this new one.
 */
void DeviceFactory::registerFactoryMethod(const DeviceCreationParams& params, DeviceFactoryMethod method)
{
	unregisterFactoryMethod(params);
	methods.insert(std::make_pair(params, method));
}

/**
 * Un-registers the device factory method referenced by the specified parameters.
 */
void DeviceFactory::unregisterFactoryMethod(const DeviceCreationParams& params)
{
	FactoryMethodsMap::iterator it = methods.find(params);
	if(it != methods.end())
	{
		methods.erase(it); // erases the element from the map
		if (methods.empty())
		{
			delete this; // A bit rock'n'roll.
		}
	}
}

/**
 * Resets the device factory: removes all factory methods.
 */
void DeviceFactory::reset()
{
	methods.clear();
}

/**
 * Creates a device supporting the specified creation parameters.
 * @throws DeviceFactory::UnsupportedDevice if no factory method matches the
 * specified creation parameters.
 */
BaseDynDevicePtr DeviceFactory::createDevice(const DeviceCreationParams& params, const uint64_t& deviceAddress)
{
	FactoryMethodsMap::const_iterator it = methods.find(params);
	if(it != methods.end())
	{
		return (*it->second)(deviceAddress); // call the device factory method that has been found.
	}
	else
	{
		throw UnsupportedDevice(params);
	}
}

//____________________________________________________________________________//

/**
 * Constructor. Initializes exception informations.
 */
DeviceFactory::UnsupportedDevice::UnsupportedDevice(const DeviceCreationParams& params)
	: params(params)
{
}

/**
 * Destructor.
 */
DeviceFactory::UnsupportedDevice::~UnsupportedDevice() throw()
{
}

void DeviceFactory::UnsupportedDevice::writeMessage(std::ostream& stream) const
{
	stream << "Unable to create a device whose parameters are: "
		"register width = " << params.registerWidth << " bytes and "
		"register type = " << params.registerType;
}

}
}
}
