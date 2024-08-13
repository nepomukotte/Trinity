/* =============================================================================
 * DeviceStoragePolicy.hpp
 * -----------------------------------------------------------------------------
 * template class mdaq::hw::DeviceStoragePolicy
 * Created on: 2 October 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef mdaq_hw_DeviceStoragePolicy_hpp_INCLUDED
#define mdaq_hw_DeviceStoragePolicy_hpp_INCLUDED

#include "AbstractControllerDevice.hpp"

namespace mdaq
{
namespace hw
{
/**
 *  This policy class delegates the storage of registers to its parent device
 *  itself.
 *  The parent device must implement the AbstractControllerDevice<T> interface.
 *  Currently, this policy is only compatible with pointer based addressing.
 */
template<typename T>
class DeviceStoragePolicy
{
public:
	typedef T* AddressType;
	typedef AbstractControllerDevice<T> DeviceType;

public:
	/// Constructor.
	/// @param device sets the parent device containing the new register.
	/// @param address sets the address of the register in the device.
	DeviceStoragePolicy(DeviceType* device = 0, AddressType address = 0);

	T read() const		{
		return device->read(regAddress);
	}
	void write(T value)	{
		device->write(regAddress, value);
	}

	AddressType address() const { return regAddress; }
	void setAddress(AddressType newAddress);

protected:
	~DeviceStoragePolicy();

private:
	DeviceStoragePolicy(const DeviceStoragePolicy&);
	DeviceStoragePolicy& operator=(const DeviceStoragePolicy&);

	DeviceType* device;		///< the device controlling the register.
	AddressType regAddress;	///< the address of the register on the device memory.
};

//____________________________________________________________________________//

/**
 * Constructor.
 */
template<typename T>
DeviceStoragePolicy<T>::DeviceStoragePolicy(DeviceType* device, AddressType address)
	: device(device), regAddress(address)
{
		device->declareRegister(regAddress);
}

/**
 * Destructor.
 */
template<typename T>
DeviceStoragePolicy<T>::~DeviceStoragePolicy()
{
		device->undeclareRegister(regAddress);
}

/**
 * Sets the address of the register.
 */
template<typename T>
void DeviceStoragePolicy<T>::setAddress(AddressType newAddress)
{
	device->undeclareRegister(regAddress);
	device->declareRegister(newAddress);
	regAddress = newAddress;
}

}
}

#endif // mdaq_hw_DeviceStoragePolicy_INCLUDED
