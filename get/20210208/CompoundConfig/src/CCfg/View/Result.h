// =============================================================================
// Result.h
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

#ifndef CCfg_View_Result_h_INCLUDED
#define CCfg_View_Result_h_INCLUDED

#include <cassert>

namespace CCfg
{
// forward declaration
class Attribute;

namespace View
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class describes the result of a view resolution.
 *
 * @author fchateau
 */
class Result
{
	Attribute* attribute;	///< the concrete attribute resulting of view resolution.
	int backtrackLevel;			///< the number of backtrackings from the theoretical attribute.

public:
	Result();
	Result(Attribute* attribute, int backtrackLevel);
	Attribute* getAttribute();
	const Attribute* getAttribute() const;
	int getBacktrackLevel() const;

	operator Attribute*();
	operator Attribute&();
	Attribute* operator ->();
	const Attribute* operator ->() const;
};

///////////////////////////////////////////////////////////////////////////////
// CCfg::CConfig inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * default constructor: creates an uninitialized view result.
 */
inline Result::Result()
	: attribute(), backtrackLevel(0)
{
}

/**
 * Constructor: creates a fully initialized view result.
 * @param attribute the concrete attribute resulting of view resolution.
 * @param backtrackLevel the number of backtracking that were necessary to find
 * the attribute.
 */
inline Result::Result(Attribute* attribute, int backtrackLevel)
	: attribute(attribute), backtrackLevel(backtrackLevel)
{
}

/**
 * Getter of the property "attribute".
 * @return Returns the concrete attribute resulting of view resolution.
 */
inline Attribute* Result::getAttribute()
{
	return attribute;
}

/**
 * Getter of the property "attribute".
 * @return Returns the concrete attribute resulting of view resolution.
 */
inline const Attribute* Result::getAttribute() const
{
	return attribute;
}

/**
 * Getter of the property "backtrackLevel".
 * @return Returns the number of backtrackings from the theoretical attribute.
 */
inline int Result::getBacktrackLevel() const
{
	return backtrackLevel;
}

/**
 * Conversion operator: returns a pointer to the concrete attribute.
 * @return Returns the concrete attribute resulting of view resolution.
 */
inline Result::operator Attribute*()
{
	return attribute;
}

/**
 * Conversion operator: returns a reference to the concrete attribute.
 * @return Returns the concrete attribute resulting of view resolution.
 */
inline Result::operator Attribute&()
{
	assert(attribute != 0);
	return *attribute;
}

/**
 * Selection operator: allows to directly access members of the concrete
 * attribute.
 * @return Returns the concrete attribute resulting of view resolution.
 */
inline Attribute* Result::operator ->()
{
	assert(attribute != 0);
	return attribute;
}

/**
 * Selection operator: allows to directly access members of the concrete
 * attribute.
 * @return Returns the concrete attribute resulting of view resolution.
 */
inline const Attribute* Result::operator ->() const
{
	assert(attribute != 0);
	return attribute;
}

}

}

#endif // CCfg_View_Result_h_INCLUDED
