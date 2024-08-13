/* =============================================================================
 * AbstractControllerDevice.hpp
 * -----------------------------------------------------------------------------
 * template class mdaq::hw::AbstractControllerDevice
 * Created on: 15 September 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef mdaq_hw_AbstractControllerDevice_h_INCLUDED
#define mdaq_hw_AbstractControllerDevice_h_INCLUDED

#include "HwControllerException.h"
#ifdef VXWORKS
#   include <../usr/h/stdint.h>
#else
#   include <stdint.h>
#endif


namespace mdaq
{
namespace hw
{
/**
 * This interface should be implemented by devices that needs to control
 * access to their registers (called "Controller Devices").
 * It is designed to be used as a base class for devices policies controlling
 * registers access.
 */
template<typename T, typename P = T*>
class AbstractControllerDevice
{
public:
	AbstractControllerDevice(uint64_t& baseAddress) : deviceAddress(baseAddress) {}
	virtual ~AbstractControllerDevice() {}

	virtual void declareRegister(P address) = 0;
	virtual void undeclareRegister(P address) = 0;
	virtual T read(P address) = 0;
	virtual void write(P address, T value) = 0;

	uint64_t& deviceAddress;

	class UndeclaredRegister;
	class RegisterAlreadyExists;
};

/**
 * This is a policy class for devices that do not need to control its registers
 * access.
 */
template<typename T, typename P = T*>
class NonControllerDevice
{
protected:
	NonControllerDevice(uint64_t& baseAddress) : deviceAddress(baseAddress) {}
	~NonControllerDevice() {}
public:
	uint64_t& deviceAddress;
};

/**
 * This exception reports that a controller device tries to access a register
 * that has not been declared.
 */
template<typename T, typename P>
class AbstractControllerDevice<T, P>::UndeclaredRegister : public HwControllerException
{
public:
	UndeclaredRegister(P address) : address(address) {}
	~UndeclaredRegister() throw() {}

	void writeMessage(std::ostream& stream) const { stream << "Register access failed: no register has been declared at address " << address; }

private:
	P address;
};

/**
 * This exception reports that a controller device cannot create a register
 * because a register already exists at the same address.
 */
template<typename T, typename P>
class AbstractControllerDevice<T, P>::RegisterAlreadyExists : public HwControllerException
{
public:
	RegisterAlreadyExists(P address) : address(address) {}
	~RegisterAlreadyExists() throw() {}

	void writeMessage(std::ostream& stream) const { stream << "Register creation failed: a register has already been declared at address " << address; }

private:
	P address;
};

}
}

#endif // mdaq_hw_AbstractControllerDevice_h_INCLUDED
