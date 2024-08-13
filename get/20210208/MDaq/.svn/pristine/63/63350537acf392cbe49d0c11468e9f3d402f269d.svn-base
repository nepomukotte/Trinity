/* =============================================================================
 * BaseDynDevice.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::hw::server::BaseDynDevice
 * class mdaq::hw::server::BaseDynDevice::RegisterNotFound
 * Created on: 23 September 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include "BaseDynDevice.h"

namespace mdaq
{
namespace hw
{
namespace server
{
/**
 * Constructor.
 */
BaseDynDevice::BaseDynDevice() : _baseAddress(0)
{

}
BaseDynDevice::BaseDynDevice(const uint64_t deviceAddress) : _baseAddress(deviceAddress)
{
}

/**
 * Destructor.
 */
BaseDynDevice::~BaseDynDevice()
{
}

/**
 * @fn void BaseDynDevice::registerList(std::vector<std::string>& list) const
 * Fills a string vector with the list of the names of the device's registers.
 */

/**
 * @fn size_t BaseDynDevice::registerCount() const
 * Returns the number of registers in the device.
 */

/**
 * @fn void BaseDynDevice::registerClear()
 * Clears the list of the device's register.
 */

/**
 * @fn BaseDynRegister&	BaseDynDevice::registerCreate(const std::string& registerName, uint64_t addressOffset)
 * Adds a new register to this device.
 * @param registerName The name of the added register.
 * @param addressOffset The offset added to the base address of the device to form the register address.
 * @return Returns a reference to the created register.
 */

/**
 * @fn BaseDynRegister* BaseDynDevice::registerSearch(const std::string& registerName)
 * Look for the register having the specified name.
 * @return The register found, 0 (null) otherwise.
 */

/**
 * @fn BaseDynRegister& BaseDynDevice::registerFind(const std::string& registerName)
 * Returns the register having the specified name.
 * @throws RegisterNotFound if no register has this name.
 */

//____________________________________________________________________________//

/**
 * @fn const std::string& BaseDynDevice::RegisterNotFound::registerName() const
 * @return The name of the missing register.
 */


/**
 * Exception constructor.
 * @param registerName The name of the missing register.
 */
BaseDynDevice::RegisterNotFound::RegisterNotFound(const std::string& registerName)
	: _registerName(registerName)
{
}

/**
 * Destructor.
 */
BaseDynDevice::RegisterNotFound::~RegisterNotFound() throw()
{
}

/**
 * Implementation of formatted message as expected by base class std::exception.
 */
void BaseDynDevice::RegisterNotFound::writeMessage(std::ostream& stream) const
{
	stream << "No register is named \"" << registerName() << '"';
}

} // namespace server
} // namespace hw
} // namespace mdaq
