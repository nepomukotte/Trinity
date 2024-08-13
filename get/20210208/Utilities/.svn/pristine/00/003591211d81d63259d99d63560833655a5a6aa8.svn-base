// =============================================================================
// BitFieldRegister.hpp
// -----------------------------------------------------------------------------
// class utl::BitFieldRegister
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef utl_BitFieldRegister_INCLUDED
#define utl_BitFieldRegister_INCLUDED

#include "utl/BitFieldHelper.hpp"
#include <boost/type_traits/remove_cv.hpp>

namespace utl
{
/**
 * This should be used as the base class of datatypes that implements a
 * bit field semantics.
 * It is designed to be used in conjunction with DECLARE_BIT and DECLARE_FIELD
 * macros. These macro must be called in the subclass body. They generate
 * accessors for each bit and fields.
 */
template<typename T>
class BitFieldRegister
{
public:
	typedef T WordType;
	typedef typename boost::remove_cv<WordType>::type WordBasicType;
	typedef utl::BitFieldHelper<WordType> RegField;

	BitFieldRegister(WordBasicType value = 0) : regValue(value) {}
	BitFieldRegister(const BitFieldRegister<T>& reg) : regValue(reg.value()) {}

	void setValue(WordBasicType val) { regValue = val; }
	WordBasicType value() const { return regValue; }

protected:
	WordType regValue;
};

#define DECLARE_RO_BIT(name, pos) \
	bool name() const { return RegField::getBit<pos>(regValue); }

#define DECLARE_WO_BIT(name, pos) \
	void set_##name(bool val) { RegField::setBit<pos>(regValue, val); }

#define DECLARE_BIT(name, pos) \
	DECLARE_RO_BIT(name, pos) \
	DECLARE_WO_BIT(name, pos)

#define DECLARE_CASTED_RO_FIELD(name, pos, width, type) \
	type name() const { return static_cast<type>(RegField::getField<pos, width>(regValue)); }

#define DECLARE_CASTED_WO_FIELD(name, pos, width, type) \
	void set_##name(type val) { RegField::setField<pos, width>(regValue, val); }

#define DECLARE_CASTED_FIELD(name, pos, width, type) \
	DECLARE_CASTED_RO_FIELD(name, pos, width, type) \
	DECLARE_CASTED_WO_FIELD(name, pos, width, type)

#define DECLARE_RO_FIELD(name, pos, width) \
	DECLARE_CASTED_RO_FIELD(name, pos, width, WordBasicType)

#define DECLARE_WO_FIELD(name, pos, width) \
	DECLARE_CASTED_WO_FIELD(name, pos, width, WordBasicType)

#define DECLARE_FIELD(name, pos, width) \
	DECLARE_RO_FIELD(name, pos, width) \
	DECLARE_WO_FIELD(name, pos, width)

}

#endif // utl_BitFieldRegister_INCLUDED
