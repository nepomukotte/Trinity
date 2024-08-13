// =============================================================================
// Range.h
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

#ifndef CCfg_Range_h_INCLUDED
#define CCfg_Range_h_INCLUDED

#include <boost/operators.hpp>
#include <string>
#include <memory>
#include <ios>

namespace CCfg
{
// forward declarations
class Value;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

class Range;
typedef std::auto_ptr<Range> RangePtr;

/**
 * This is an abstract class defining the concept of range: a restriction on
 * the values a variable can have.
 *
 * @author fchateau
 */
class Range : boost::equality_comparable<Range>
{
public:
	Range();
	virtual ~Range();

	virtual bool contains(const Value& value) const = 0;
	virtual bool equals(const Range& rvalue) const = 0;
	virtual RangePtr clone() const = 0;
	virtual std::string toString() const = 0;
	virtual void fromString(const std::string& str) = 0;
	virtual void write(std::ostream& stream) const = 0;

	bool operator==(const Range& rvalue) const;

	static RangePtr createFromString(const std::string& str);
	static bool equivalent(const Range* lvalue, const Range* rvalue);
};

std::ostream& operator <<(std::ostream& stream, const Range& range);



////////////////////////////////////////////////////////////////////////////////
// CCfg::Range inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Equality comparison operator.
 */
inline bool Range::operator==(const Range& rvalue) const
{
	return equals(rvalue);
}

/**
 * Stream insertion operator.
 */
inline std::ostream& operator <<(std::ostream& stream, const Range& range)
{
	range.write(stream);
	return stream;
}

}

#endif // CCfg_Range_h_INCLUDED
