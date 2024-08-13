// =============================================================================
// Range.cpp
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

#include "Range.h"
#include "Grammar/RangeGrammar.h"

namespace CCfg
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Range class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
Range::Range()
{
}

/**
 * Destructor.
 */
Range::~Range()
{
}

/**
 * @fn bool Range::contains(const Value& value) const
 * Returns whether the range contains the specified value.
 */

/**
 * @fn bool Range::equals(const Range& rvalue) const
 * Returns whether two ranges are equals.
 */

/**
 * @fn RangePtr Range::clone() const
 * Allocates a copy of this range.
 */

/**
 * @fn std::string Range::toString() const
 * Converts the range to a string format.
 */

/**
 * @fn void CCfg::Range::fromString(const std::string& str)
 * Parses the string representation of a range into *this.
 */

/**
 * @fn void Range::write(std::ostream& stream) const
 * Writes the string representation of this range on a stream.
 */

/**
 * Creates a range from a string compatible with a range format.
 */
RangePtr Range::createFromString(const std::string& str)
{
	return Grammar::RangeGrammar::fromString(str);
}

/**
 * Compares 2 range pointers for equality.
 * They are equivalent if they are equal or if the ranges they reference are the same.
 */
bool Range::equivalent(const Range* lvalue, const Range* rvalue)
{
	return (lvalue == rvalue) or
		(lvalue != 0 and rvalue != 0 and *lvalue == *rvalue);
}

}
