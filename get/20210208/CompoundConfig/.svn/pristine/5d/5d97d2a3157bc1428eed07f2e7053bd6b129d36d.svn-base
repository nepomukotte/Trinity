// =============================================================================
// Interval.h
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

#ifndef CCfg_Interval_h_INCLUDED
#define CCfg_Interval_h_INCLUDED

#include "Range.h"
#include "Value.h"

namespace CCfg
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class defines an interval: a range which enclose numbers between two boundaries.
 * It can define either integer, hexadecimal, real or string (lexicographical
 * order) intervals.
 *
 * @author fchateau
 */
class Interval : public Range
{
public:
	// enumerations
	enum BoundMode {IncludeBound, ///< the value of the bound is included in the interval (open interval).
		ExcludeBound,	///< the value of the bound is included for the interval (closed interval).
		Unbounded ///< it represents infinity: the value of the bound is ignored and tests are not performed against this bound.
	};

	// static constants
	static const char LBOUND_INC = '[';
	static const char LBOUND_EXC = ']';
	static const char UBOUND_INC = ']';
	static const char UBOUND_EXC = '[';
	static const char SEPARATOR = ',';
	static const char RESERVED_CHARACTERS[]; // concatenation of previous characters
	static const char MINUS_INFINITY[];
	static const char PLUS_INFINITY[];

	Interval();
	Interval(const Value& lower, const Value& upper, BoundMode lowerMode = IncludeBound, BoundMode upperMode = IncludeBound);
	Interval(Int lower, Int upper, BoundMode lowerMode = IncludeBound, BoundMode upperMode = IncludeBound);
	Interval(Hex lower, Hex upper, BoundMode lowerMode = IncludeBound, BoundMode upperMode = IncludeBound);
	Interval(Real lower, Real upper, BoundMode lowerMode = IncludeBound, BoundMode upperMode = IncludeBound);
	Interval(const String& lower, const String& upper, BoundMode lowerMode = IncludeBound, BoundMode upperMode = IncludeBound);
	Interval(const Interval& rvalue);
	virtual ~Interval();
	Interval& operator=(const Interval& rvalue);

	virtual bool contains(const Value& value) const;
	virtual bool equals(const Range& rvalue) const;
	virtual RangePtr clone() const;
	virtual std::string toString() const;
	virtual void fromString(const std::string& str);
	virtual void write(std::ostream& stream) const;

	const Value& getLowerBound() const;
	const Value& getUpperBound() const;
	void setLowerBound(const Value& lower);
	void setUpperBound(const Value& upper);

	BoundMode getLowerBoundMode() const;
	BoundMode getUpperBoundMode() const;
	void setLowerBoundMode(BoundMode mode);
	void setUpperBoundMode(BoundMode mode);

private:
	Value lowerBound;	///< the inferior boundary of the interval.
	Value upperBound;	///< the superior boundary of the interval.
	BoundMode lowerBoundMode;	///< the mode of the lower bound: inclusive, exclusive, unbounded.
	BoundMode upperBoundMode;	///< the mode of the upper bound
};

std::ostream& operator <<(std::ostream& stream, const Interval& interval);

////////////////////////////////////////////////////////////////////////////////
// CCfg::Interval inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Returns the inferior boundary.
 */
inline const Value& Interval::getLowerBound() const
{
	return lowerBound;
}

/**
 * Returns the superior boundary.
 */
inline const Value& Interval::getUpperBound() const
{
	return upperBound;
}

/**
 * Redefines the inferior boundary.
 */
inline void Interval::setLowerBound(const Value& lower)
{
	lowerBound = lower;
}

/**
 * Redefines the superior boundary.
 */
inline void Interval::setUpperBound(const Value& upper)
{
	upperBound = upper;
}


/**
 * Returns the mode of the lower bound.
 */
inline Interval::BoundMode Interval::getLowerBoundMode() const
{
	return lowerBoundMode;
}

/**
 * Returns the mode of the upper bound
 */
inline Interval::BoundMode Interval::getUpperBoundMode() const
{
	return upperBoundMode;
}

/**
 * Defines the mode of the lower bound.
 */
inline void Interval::setLowerBoundMode(BoundMode mode)
{
	lowerBoundMode = mode;
}

/**
 * Defines the mode of the upper bound.
 */
inline void Interval::setUpperBoundMode(BoundMode mode)
{
	upperBoundMode = mode;
}

/**
 * Stream insertion operator.
 */
inline std::ostream& operator <<(std::ostream& stream, const Interval& interval)
{
	interval.write(stream);
	return stream;
}



}

#endif // CCfg_Interval_h_INCLUDED
