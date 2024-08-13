/**
 * @file EnumStringMapper.hpp
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

#ifndef utl_EnumStringMapper_hpp_INCLUDED
#define utl_EnumStringMapper_hpp_INCLUDED

#include <utl/Exception.h>
#include <utl/compat/ArrayHelper.hpp>
#include <sstream>
#include <typeinfo>

namespace utl
{

/**
 * This class allows conversions between enumerated values and their string representation.
 * If the enum values are not contiguous and starting from zero, the list of values must
 * be declared once in a .cpp using the ENUM_LOOKUP_DEF macro.
 * In both case, the ENUM_STRINGS_DEF macro must be used to declare the string representations
 * of each value.
 */
template<class Enum>
struct EnumStringMapper
{
	static Enum fromString(const std::string& str);
	static const char* toString(Enum e);
};

template<class Enum>
struct EnumIndexMapper
{
	static size_t indexOf(Enum value);
	static Enum at(size_t index);

private:
	struct Impl;
};

class EnumStringNotFound : public Exception
{
public:
	EnumStringNotFound(const std::type_info& eti, const std::string& valueStr);
	~EnumStringNotFound() throw();

	void writeMessage(std::ostream& stream) const;

private:
	std::string valueStr;
	std::string enumName;
};

class EnumMapperBadDefinition : public Exception
{
public:
	explicit EnumMapperBadDefinition(const std::type_info& eti);
	EnumMapperBadDefinition(const EnumMapperBadDefinition& rvalue);
	~EnumMapperBadDefinition() throw();
	EnumMapperBadDefinition& operator=(const EnumMapperBadDefinition& rvalue);

	template<class T> EnumMapperBadDefinition& operator <<(T data)
	{
		msgBuilder << data;
		return *this;
	}

	void writeMessage(std::ostream& stream) const;

private:
	std::ostringstream msgBuilder;
};

} /* namespace utl */

#endif /* utl_EnumStringMapper_hpp_INCLUDED */
