/**
 * @file EnumStringMapper.ipp
 * @date 21 juil. 2014
 * @author fchateau
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Frederic Chateau
 *
 * This file is part of the Utilities software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
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
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#ifndef utl_EnumStringMapper_ipp_INCLUDED
#define utl_EnumStringMapper_ipp_INCLUDED

#include "EnumStringMapper.hpp"

#define ENUM_STRINGS_DEF(enumtype) \
	template struct ::utl::EnumStringMapper< enumtype >; \
	template struct ::utl::EnumIndexMapper< enumtype >; \
	template<> const char* const ::utl::EnumStrings< enumtype >::values[]

#define ENUM_LOOKUP_DEF(enumtype) \
	namespace utl { template<> struct EnumIndexMapper< enumtype >::Impl : EnumLUTMapping< enumtype > {}; } \
	template<> enumtype const ::utl::EnumLUTMapping< enumtype >::lookUpTable[]

namespace utl
{

template<class E>
struct EnumIdentityMapping
{
	static size_t indexOf(E value) { return value; }
	static E at(size_t index) { return static_cast<E>(index); }
};

template<class E>
struct EnumLUTMapping
{
	static size_t indexOf(E value)
	{
		size_t index = array::find_index(array::begin(lookUpTable), array::end(lookUpTable), value);
		if(index >= array::size(lookUpTable))
		{
			throw EnumMapperBadDefinition(typeid(E)) << "cannot find the value " << value << " in the lookup table";
		}
		return index;
	}
	static E at(size_t index)
	{
		if(index >= array::size(lookUpTable))
		{
			throw EnumMapperBadDefinition(typeid(E)) << "no element #" << index << " in the lookup table";
		}
		return lookUpTable[index];
	}
private:
	static const E lookUpTable[];
};

template<class E>
struct EnumStrings
{
	static const char* const values[];
	static const char* const * begin() { return array::begin(values); }
	static const char* const * end()   { return array::end  (values); }
	static size_t size()               { return array::size(values); }
};

template<class E>
struct EnumIndexMapper<E>::Impl : EnumIdentityMapping<E>
{
};

// =============================================================================
// utl::EnumIndexMapper inline methods body
// -----------------------------------------------------------------------------

template<class Enum>
size_t EnumIndexMapper<Enum>::indexOf(Enum value)
{
	return Impl::indexOf(value);
}

template<class Enum>
Enum EnumIndexMapper<Enum>::at(size_t index)
{
	return Impl::at(index);
}

// =============================================================================
// utl::EnumStringMapper inline methods body
// -----------------------------------------------------------------------------

template<class Enum>
Enum EnumStringMapper<Enum>::fromString(const std::string& str)
{
	size_t index = array::find_index(EnumStrings<Enum>::begin(), EnumStrings<Enum>::end(), str);
	if(index >= EnumStrings<Enum>::size())
	{
		throw EnumStringNotFound(typeid(Enum), str);
	}
	return EnumIndexMapper<Enum>::at(index);
}

template<class Enum>
const char* EnumStringMapper<Enum>::toString(Enum e)
{
	size_t index = EnumIndexMapper<Enum>::indexOf(e);
	if(index >= EnumStrings<Enum>::size())
	{
		throw EnumMapperBadDefinition(typeid(Enum)) << "cannot find the name of the enumerated value " << e;
	}
	return EnumStrings<Enum>::values[index];
}

} /* namespace utl */

#endif /* utl_EnumStringMapper_ipp_INCLUDED */
