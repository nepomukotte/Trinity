/* =============================================================================
 * BaseDynDevice.h
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

#ifndef mdaq_hw_server_BaseDynDevice_h_INCLUDED
#define mdaq_hw_server_BaseDynDevice_h_INCLUDED

#include "BaseDynRegister.h"
#include "../HwControllerException.h"
#include <vector>
#include <string>
#include <memory>
#include <iostream>

namespace mdaq
{
namespace hw
{
namespace server
{
/**
 * This is the base class for accessing and configuring a device dynamically.
 */
class BaseDynDevice
{
public:
	class RegisterNotFound;

	BaseDynDevice();
	BaseDynDevice(const uint64_t deviceAddress);
	virtual ~BaseDynDevice();

	void setBaseAddress(uint64_t address) {_baseAddress = address;}
	uint64_t baseAddress() const { return _baseAddress; }

	virtual void	registerList(std::vector<std::string>& list) const = 0;
	virtual size_t	registerCount() const = 0;
	virtual void	registerClear() = 0;

	virtual BaseDynRegister& registerCreate	(const std::string& registerName, uint64_t addressOffset) = 0;
	virtual BaseDynRegister* registerSearch	(const std::string& registerName) = 0;
	virtual BaseDynRegister& registerFind	(const std::string& registerName) = 0;
	uint64_t _baseAddress;
};

typedef std::auto_ptr<BaseDynDevice> BaseDynDevicePtr;

/**
 * Exception thrown by classes deriving from BaseDynDevice when no register
 * exists for the specified name.
 */
class BaseDynDevice::RegisterNotFound : public HwControllerException
{
public:
	RegisterNotFound(const std::string& registerName);
	~RegisterNotFound() throw();
	const std::string& registerName() const { return _registerName; }

private:
	void writeMessage(std::ostream& stream) const;
	std::string _registerName; ///< the name of the missing register.
};

} // namespace server
} // namespace hw
} // namespace mdaq

#endif // mdaq_hw_server_BaseDynDevice_h_INCLUDED
