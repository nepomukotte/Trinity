// =============================================================================
// TypeInfo.cpp
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

#include "TypeInfo.h"
#include <string.h>

namespace CCfg
{

struct TypeData
{
	AttributeTypeId attrId;
	const char* const name;
	const char* const defValue;
} const types[VALUETYPE_COUNT] =
{
	{INT_SCALAR, "integer", "0"},
	{HEX_SCALAR, "hexadecimal", "0x0"},
	{REAL_SCALAR, "real number", "0.0"},
	{BOOL_SCALAR, "boolean", "false"},
	{STRING_SCALAR, "string", ""}
};

////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::TypeInfo class methods body
////////////////////////////////////////////////////////////////////////////////


/**
 * Returns the id of a type given its display name.
 */
ValueTypeId TypeInfo::lookup(const char* name)
{
	for(unsigned i = 0; i < sizeof types / sizeof types[0]; ++i)
	{
		if(strcmp(name, types[i].name)	== 0)
		{
			return static_cast<ValueTypeId>(i);
		}
	}
	return static_cast<ValueTypeId>(VALUETYPE_COUNT);
}

/**
 * Returns the id of a type (ValueTypeId) given its AttributeTypeId.
 */
ValueTypeId TypeInfo::lookup(AttributeTypeId attrId)
{
	if(attrId == OBJECT)
	{
		return VALUETYPE_INVALID;
	}
	else
	{
		return static_cast<ValueTypeId>(attrId - 1);
	}
}


/**
 * Returns the AttributeTypeId of type given its id.
 */
AttributeTypeId TypeInfo::getAttributeTypeId(ValueTypeId id)
{
	return types[id].attrId;
}

/**
 * Returns the display name of a type given its id.
 */
const char* TypeInfo::getDisplayName(ValueTypeId id)
{
	return types[id].name;
}

const char* TypeInfo::getDisplayName(AttributeTypeId id)
{
	if(id == OBJECT)
	{
		return "object";
	}
	else
	{
		return types[id-1].name;
	}
}

const char* TypeInfo::getDefaultValue(ValueTypeId id)
{
	return types[id].defValue;
}

}
