/* =============================================================================
 * DynDevice.hpp
 * -----------------------------------------------------------------------------
 * template class mdaq::hw::server::DynDevice
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

#ifndef mdaq_hw_DynDevice_hpp_INCLUDED
#define mdaq_hw_DynDevice_hpp_INCLUDED

#include "BaseDynDevice.h"
#include "DynRegister.hpp"
#include "../MemStoragePolicy.hpp"
#include "../DeviceStoragePolicy.hpp"
#include "../AbstractControllerDevice.hpp"
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/type_traits/is_member_function_pointer.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

namespace mdaq
{
namespace hw
{
namespace server
{
/**
 * The class represents a hardware device as set of registers.
 */
template<typename T,
	template <class Data> class S = MemStoragePolicy,
	template <class Data, class Pointer> class C = NonControllerDevice>
class DynDevice : public BaseDynDevice, public C<T, typename S<T>::AddressType>
{
public:
	/// The word type of the registers of this device
	typedef T WordType;
	/// How the register is stored and accessed.
	typedef S<T> StoragePolicy;
	/// The type of addresses.
	typedef typename StoragePolicy::AddressType AddressType;
	/// How the device controls the storage of registers.
	typedef C<T, AddressType> ControllerPolicy;
	/// The register type associated to this device type
	typedef DynRegister<T, S> Register;

	/// Constructor
	DynDevice(uint64_t deviceAddress) : BaseDynDevice(deviceAddress), ControllerPolicy(BaseDynDevice::_baseAddress) {}

	void		registerList(std::vector<std::string>& list) const;
	size_t		registerCount() const	{ return regs.size(); }
	void		registerClear()			{ regs.clear(); }

	Register&	registerCreate	(const std::string& registerName, uint64_t addressOffset);
	Register*	registerSearch	(const std::string& registerName);
	Register&	registerFind	(const std::string& registerName);


	/// Easy register access operator.
	Register& operator[] (const std::string& registerName) { return registerFind(registerName); }

	static BaseDynDevicePtr factoryMethod(const uint64_t& deviceAddress)
	{
		return BaseDynDevicePtr(new DynDevice<T,S,C> (deviceAddress));
	}


protected:
	/// All the registers of this device mapped by their name.
	typedef boost::ptr_map<std::string, Register> RegisterMap;

	RegisterMap regs; ///< Registers of this device, mapped by their name.

	template<typename R>
	Register& registerCreateImpl(const std::string& registerName,
		typename boost::enable_if<boost::is_base_of<DeviceStoragePolicy<T>, R>, uint64_t>::type addressOffset)
	{
		Register* newReg = new Register(registerName, this, reinterpret_cast<AddressType>(addressOffset));
		return *regs.insert(const_cast<std::string&>(registerName), newReg).first->second;
	}

	template<typename R>
	Register& registerCreateImpl(const std::string& registerName,
		typename boost::enable_if<boost::is_base_of<MemStoragePolicy<T>, R>, uint64_t>::type addressOffset)
	{
		Register* newReg = new Register(registerName, this, reinterpret_cast<AddressType>(addressOffset));
		return *regs.insert(const_cast<std::string&>(registerName), newReg).first->second;
	}

	template<typename R>
	Register& registerCreateImpl(const std::string& registerName,
		typename boost::enable_if<boost::is_base_of<NestedStoragePolicy<T>, R>, uint64_t>::type addressOffset)
	{
		Register* newReg = new Register(registerName);
		newReg->setStorageAddress(addressOffset);
		return *regs.insert(const_cast<std::string&>(registerName), newReg).first->second;
	}
};

//____________________________________________________________________________//

/*
 * (see inherited doc)
 */
template<typename T, template <class Data> class S, template <class Data, class Pointer> class C>
void
DynDevice<T, S, C>::registerList(std::vector<std::string>& list) const
{
	list.clear();
	for(typename RegisterMap::const_iterator it = regs.begin(), end = regs.end();
		it != end; ++it)
	{
		list.push_back(it->first);
	}
}

/*
 * (see inherited doc)
 */
template<typename T, template <class Data> class S, template <class Data, class Pointer> class C>
typename DynDevice<T, S, C>::Register&
DynDevice<T, S, C>::registerCreate(const std::string& registerName, uint64_t addressOffset)
{
	return registerCreateImpl<Register>(registerName, addressOffset);
}

/*
 * (see inherited doc)
 */
template<typename T, template <class Data> class S, template <class Data, class Pointer> class C>
typename DynDevice<T, S, C>::Register*
DynDevice<T, S, C>::registerSearch(const std::string& registerName)
{
	typename RegisterMap::iterator it = regs.find(registerName);
	return it == regs.end() ? 0 : it->second;
}

/*
 * (see inherited doc)
 */
template<typename T, template <class Data> class S, template <class Data, class Pointer> class C>
typename DynDevice<T, S, C>::Register&
DynDevice<T, S, C>::registerFind(const std::string& registerName)
{
	typename RegisterMap::iterator it = regs.find(registerName);
	if(it == regs.end())
	{
		throw RegisterNotFound(registerName);
	}
	return *it->second;
}

} // namespace server
} // namespace hw
} // namespace mdaq

#endif // mdaq_hw_DynDevice_hpp_INCLUDED
