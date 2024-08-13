/* =============================================================================
 * DynRegister.hpp
 * -----------------------------------------------------------------------------
 * template class mdaq::hw::server::DynRegister
 * template class mdaq::hw::server::DynRegister::FieldOutOfRange
 * Created on: 28 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef mdaq_hw_server_DynRegister_hpp_INCLUDED
#define mdaq_hw_server_DynRegister_hpp_INCLUDED

#include "../SmartRegister.hpp"
#include "../MemStoragePolicy.hpp"
//#include "../NestedStoragePolicy.hpp"
#include "BaseDynRegister.h"
#include <utl/BitFieldHelper.hpp>
#include <limits>

#undef min
#undef max

namespace mdaq
{
namespace hw
{
namespace server
{
/**
 * This class represents a hardware register optionally composed of bit-fields.
 * Bit-fields should be given a name to reference them.
 */
template<typename T, template <class Type> class S = MemStoragePolicy>
class DynRegister : public SmartRegister<T, S, NonObservableValue, BaseDynRegister>
{
public:
	typedef SmartRegister<T, S, NonObservableValue, BaseDynRegister> Register;
	typedef typename Register::AddressType AddressType;
	typedef typename Register::DeviceType DeviceType;

	DynRegister(const std::string& name) : Register(), name(name) {}
	DynRegister(const std::string& name, AddressType address) : Register(address), name(name) {}
	DynRegister(const std::string& name, DeviceType* device, AddressType address = 0) : Register(device, address), name(name) {}
	~DynRegister() {}

	void setStorageAddress(uint64_t address);
	uint64_t storageAddress() const;

	T		fieldRead(const std::string& name) const;
	void	fieldWrite(const std::string& name, T value);

	void	setValueInt64(int64_t value);
	int64_t getValueInt64() const;

	int64_t	fieldReadInt64(const std::string& name) const;
	void	fieldWriteInt64(const std::string& name, int64_t value);

	class OutOfRange;

private:
	std::string name;
};

/**
 * This exception is thrown whenever the program tries to set a bit-field with
 * a value over its size.
 */
template<typename T, template <class Type> class S>
class DynRegister<T, S>::OutOfRange : public BaseDynRegister::OutOfRange
{
public:
	OutOfRange(const std::string& itemName, BaseDynRegister::OutOfRangeType type, T value, T minValue, T maxValue);
	~OutOfRange() throw() {}

private:
	void writeMessage(std::ostream& stream) const;

	T value;
	T minValue;
	T maxValue;
};

//____________________________________________________________________________//

template<typename T, template <class Type> class S>
inline void DynRegister<T, S>::setStorageAddress(uint64_t address)
{
	Register::StoragePolicy::setAddress(reinterpret_cast<typename DynRegister::AddressType>(address));
}

template<typename T, template <class Type> class S>
inline uint64_t DynRegister<T, S>::storageAddress() const
{
	return reinterpret_cast< uint64_t >(Register::StoragePolicy::address());
}

/**
 * Reads the value of the specified bit-field.
 */
template<typename T, template <class Type> class S>
inline T DynRegister<T, S>::fieldRead(const std::string& name) const
{
	typename BaseDynRegister::BitFieldsMap::const_iterator it = this->fields.find(name);
	if(it != this->fields.end())
	{
		const BaseDynRegister::BitFieldParameters& field = it->second;
		return utl::BitFieldHelper<T>::getField(this->value(), field.offset, field.width);
	}
	else
	{
		throw BaseDynRegister::FieldNotFound(name);
	}
}

/**
 * Writes the value of the specified bit-field.
 */
template<typename T, template <class Type> class S>
inline void DynRegister<T, S>::fieldWrite(const std::string& name, T value)
{
	typename BaseDynRegister::BitFieldsMap::const_iterator it = this->fields.find(name);
	if(it != this->fields.end())
	{
		const BaseDynRegister::BitFieldParameters& field = it->second;
		T maxValue = (1 << field.width) - 1;
		if (value > maxValue)
		{
			throw OutOfRange(name, BaseDynRegister::FieldOutOfRange, value, 0, maxValue);
		}
		T regValue = this->value();
		utl::BitFieldHelper<T>::setField(regValue, field.offset, field.width, value);
		this->setValue(regValue);
	}
	else
	{
		throw BaseDynRegister::FieldNotFound(name);
	}
}

/*
 * (see inherited doc)
 */
template<typename T, template <class Type> class S>
inline void DynRegister<T, S>::setValueInt64(int64_t value)
{
	int64_t minValue = std::numeric_limits<T>::min();
	int64_t maxValue = std::numeric_limits<T>::max();
	if( static_cast<uint64_t>(value) > static_cast<uint64_t>(maxValue)
			or static_cast<uint64_t>(value) < static_cast<uint64_t>(minValue) )
	{
		throw OutOfRange(name, BaseDynRegister::RegisterOutOfRange, value, minValue, maxValue);
	}
	this->setValue(static_cast<T>(value));
}

/*
 * (see inherited doc)
 */
template<typename T, template <class Type> class S>
inline int64_t DynRegister<T, S>::getValueInt64() const
{
	return static_cast<int64_t>(this->value());
}

/*
 * (see inherited doc)
 */
template<typename T, template <class Type> class S>
int64_t DynRegister<T, S>::fieldReadInt64(const std::string& name) const
{
	return static_cast<int64_t>(fieldRead(name));
}

/*
 * (see inherited doc)
 */
template<typename T, template <class Type> class S>
void DynRegister<T, S>::fieldWriteInt64(const std::string& name, int64_t value)
{
	return this->fieldWrite(name, static_cast<T>(value));
}

//____________________________________________________________________________//

/**
 * Exception constructor.
 */
template<typename T, template <class Type> class S>
inline DynRegister<T, S>::OutOfRange::OutOfRange(const std::string& itemName, BaseDynRegister::OutOfRangeType type, T value, T minValue, T maxValue)
	: BaseDynRegister::OutOfRange(itemName, type),
	value(value),
	minValue(minValue),
	maxValue(maxValue)
{
}

/**
 * Implementation of formatted message as expected by base class std::exception.
 */
template<typename T, template <class Type> class S>
inline void DynRegister<T, S>::OutOfRange::writeMessage(std::ostream& stream) const
{
	stream << "Error setting the ";
	if(type == BaseDynRegister::RegisterOutOfRange)
	{
		stream << "register";
	}
	else
	{
		stream << "field";
	}
	stream <<" \"" << itemName << "\" value to " << value
		<< " (valid range is " << minValue << " to " << maxValue << ')';
}

}
}
}

#endif // mdaq_hw_server_DynRegister_hpp_INCLUDED
