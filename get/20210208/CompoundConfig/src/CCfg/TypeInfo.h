// =============================================================================
// TypeInfo.h
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef CCfg_TypeInfo_h_INCLUDED
#define CCfg_TypeInfo_h_INCLUDED

#include <string>

namespace CCfg
{

// enums

enum ValueTypeId
{
	INTEGER,
	HEXADECIMAL,
	REAL,
	BOOLEAN,
	STRING,
	VALUETYPE_INVALID
};

#define VALUETYPE_COUNT VALUETYPE_INVALID

enum AttributeTypeId
{
	OBJECT,
	INT_SCALAR,
	HEX_SCALAR,
	REAL_SCALAR,
	BOOL_SCALAR,
	STRING_SCALAR
};

// typedefs
typedef int			Int;
typedef unsigned	Hex;
typedef double		Real;
typedef bool		Bool;
typedef std::string String;


////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class provides informations about types used in scalar attributes.
 *
 * @author fchateau
 */
class TypeInfo
{
	TypeInfo();
public:
	static ValueTypeId lookup(const char* name);
	static ValueTypeId lookup(AttributeTypeId attrId);

	static AttributeTypeId getAttributeTypeId(ValueTypeId id);
	static const char* getDisplayName(ValueTypeId id);
	static const char* getDisplayName(AttributeTypeId id);
	static const char* getDefaultValue(ValueTypeId id);
};

/**
 * This class realizes a compile time mapping between a type and its enum id.
 * It is mainly useful for template functions.
 * example of use: ValueTypeId id = Id<T>::val;
 *
 * @author fchateau
 */
template<typename T>
class Id
{
};

#define DECL_TYPEMAP(type, value)\
	template<>\
	struct Id<type>\
	{\
		static const ValueTypeId val = value;\
	};

DECL_TYPEMAP(Int, INTEGER)
DECL_TYPEMAP(Hex, HEXADECIMAL)
DECL_TYPEMAP(Real, REAL)
DECL_TYPEMAP(Bool, BOOLEAN)
DECL_TYPEMAP(std::string, STRING)

}

#endif // CCfg_TypeInfo_h_INCLUDED
