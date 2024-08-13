/**
 * @file InplaceParser.ipp
 * @date 5 juin 2013
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

#ifndef utl_io_InplaceParser_ipp_INCLUDED
#define utl_io_InplaceParser_ipp_INCLUDED

#include "InplaceParser.hpp"
#include <boost/utility/enable_if.hpp>

namespace utl
{
namespace io
{
// =============================================================================
// utl::io::InplaceParser::has_checker class definition
// -----------------------------------------------------------------------------

template<class T>
struct InplaceParser::has_checker
{
	typedef char Yes[2];
	typedef char No[1];

	template<class C, C> class Test;
	template<class U> static Yes& check(Test<bool(*)(const InplaceParser&, const U*), &isDataField>*);
	template<class U> static No& check(...);

	const static bool value = sizeof(check<T>(0)) == sizeof(Yes);
};

// =============================================================================
// utl::io::InplaceParser class methods body
// -----------------------------------------------------------------------------

inline InplaceParser::InplaceParser()
	: bufCurPtr(0), bufEndPtr(0)
{
}

inline InplaceParser::InplaceParser(const ConstBuffer& buf)
	: bufCurPtr(reinterpret_cast<const char*>(buf.data)),
	  bufEndPtr(bufCurPtr + buf.maxBytes)
{
}

inline InplaceParser::InplaceParser(const InplaceParser& p)
	: bufCurPtr(p.bufCurPtr),
	  bufEndPtr(p.bufEndPtr)
{
}

inline InplaceParser& InplaceParser::operator=(const InplaceParser& val)
{
	bufCurPtr = val.bufCurPtr;
	bufEndPtr = val.bufEndPtr;
	return *this;
}

template<class T>
bool InplaceParser::bind(const T*& val, size_t n)
{
	if(bufEndPtr >= bufCurPtr + sizeof(T) * n)
	{
		val = reinterpret_cast<const T*>(bufCurPtr);
		bufCurPtr += sizeof(T) * n;
		return true;
	}
	else
	{
		return false;
	}
}

template<class T>
bool InplaceParser::tryBind(const T*& val, size_t n) const
{
	if(bufEndPtr >= bufCurPtr + sizeof(T) * n)
	{
		val = reinterpret_cast<const T*>(bufCurPtr);
		return true;
	}
	else
	{
		return false;
	}
}

template<class T>
bool InplaceParser::ignore(size_t n)
{
	if(bufEndPtr >= bufCurPtr + sizeof(T) * n)
	{
		bufCurPtr += sizeof(T) * n;
		return true;
	}
	else
	{
		return false;
	}
}

inline bool InplaceParser::commit(InplaceParser& dest) const
{
	dest = *this;
	return true;
}

inline bool InplaceParser::isEndOfBuffer() const
{
	return cur() == end();
}

template<class DataField>
bool parse_1_impl(
		typename boost::enable_if_c<
			InplaceParser::has_checker<DataField>::value,
			InplaceParser&
		>::type buf,
		ValueRef<DataField>& dest)
{
	return isDataField(buf, dest.get()) && buf.bind(dest.get());
}

template<class DataField>
bool parse_1_impl(
		typename boost::disable_if_c<
			InplaceParser::has_checker<DataField>::value,
			InplaceParser&
		>::type buf,
		ValueRef<DataField>& dest)
{
	return buf.bind(dest.get());
}

template<class DataField>
bool InplaceParser::parse_1(ValueRef<DataField>& dest)
{
	return parse_1_impl(*this, dest);
}

/**
 * Binds N elements without checking data prefix.
 */
template<class DataField>
bool InplaceParser::parse_N(ArrayRef<DataField>& dest, unsigned count)
{
	if(bind(dest.get(), count))
	{
		dest.count() = count;
		return true;
	}
	else
	{
		return false;
	}
}

template<class DataField>
bool InplaceParser::parse_0_1(ValueRef<DataField>& dest)
{
	parse_1(dest);
	return true;
}

template<class DataField>
bool InplaceParser::parse_0_N(ArrayRef<DataField>& dest)
{
	parse_1_N(dest);
	return true;
}

template<class DataField>
bool InplaceParser::parse_1_N(ArrayRef<DataField>& dest)
{
	if(isDataField(*this, dest.get()))
	{
		bind(dest.get());
		++dest.count();

		while(isDataField(*this, dest.get()))
		{
			const DataField* nextField;
			bind(nextField);
			++dest.count();
		}
		return true;
	}
	else
	{
		return false;
	}
}

template<class BindingField>
bool InplaceParser::parse_0_N(std::vector<BindingField>& dest)
{
	parse_1_N(dest);
	return true;
}

template<class BindingField>
bool InplaceParser::parse_1_N(std::vector<BindingField>& dest)
{
	while(cur() < end())
	{
		BindingField element;
		if(parse_1(element))
		{
			dest.push_back(element);
		}
		else
		{
			break;
		}
	}
	return dest.size() > 0;
}

} /* namespace io */
} /* namespace utl */

#endif /* utl_io_InplaceParser_ipp_INCLUDED */
