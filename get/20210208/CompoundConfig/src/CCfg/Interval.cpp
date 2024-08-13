// =============================================================================
// Interval.cpp
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

#include "Interval.h"
#include "Grammar/IntervalGrammar.h"
#include <sstream>

namespace CCfg
{

////////////////////////////////////////////////////////////////////////////////
// CCfg::Interval class methods body
////////////////////////////////////////////////////////////////////////////////

const char Interval::LBOUND_INC;
const char Interval::LBOUND_EXC;
const char Interval::UBOUND_INC;
const char Interval::UBOUND_EXC;
const char Interval::SEPARATOR;
const char Interval::RESERVED_CHARACTERS[] = ",[]";
const char Interval::MINUS_INFINITY[] = "-inf";
const char Interval::PLUS_INFINITY[] = "+inf";

/**
 * Default constructor.
 */
Interval::Interval()
	: lowerBound(), upperBound(),
	lowerBoundMode(IncludeBound), upperBoundMode(IncludeBound)
{
}

/**
 * Full constructor with generic values.
 */
Interval::Interval(const Value& lower, const Value& upper,
	BoundMode lowerBoundMode, BoundMode upperBoundMode)
	: lowerBound(lower), upperBound(upper),
	lowerBoundMode(lowerBoundMode), upperBoundMode(upperBoundMode)
{
}

/**
 * Full constructor with integer values.
 */
Interval::Interval(Int lower, Int upper,
	BoundMode lowerBoundMode, BoundMode upperBoundMode)
	: lowerBound(lower), upperBound(upper),
	lowerBoundMode(lowerBoundMode), upperBoundMode(upperBoundMode)
{
}

/**
 * Full constructor with hexadecimal values.
 */
Interval::Interval(Hex lower, Hex upper,
	BoundMode lowerBoundMode, BoundMode upperBoundMode)
	: lowerBound(lower), upperBound(upper),
	lowerBoundMode(lowerBoundMode), upperBoundMode(upperBoundMode)
{
}

/**
 * Full constructor with real values.
 */
Interval::Interval(Real lower, Real upper,
	BoundMode lowerBoundMode, BoundMode upperBoundMode)
	: lowerBound(lower), upperBound(upper),
	lowerBoundMode(lowerBoundMode), upperBoundMode(upperBoundMode)
{
}

/**
 * Full constructor with string values.
 */
Interval::Interval(const String& lower, const String& upper,
	BoundMode lowerBoundMode, BoundMode upperBoundMode)
	: lowerBound(lower), upperBound(upper),
	lowerBoundMode(lowerBoundMode), upperBoundMode(upperBoundMode)
{
}

/**
 * Copy constructor.
 */
Interval::Interval(const Interval& rvalue)
	: Range(rvalue), lowerBound(rvalue.lowerBound), upperBound(rvalue.upperBound),
	lowerBoundMode(rvalue.lowerBoundMode), upperBoundMode(rvalue.upperBoundMode)
{
}

/**
 * Destructor.
 */
Interval::~Interval()
{
}

/**
 * Assignment operator.
 */
Interval& Interval::operator=(const Interval& rvalue)
{
	lowerBound = rvalue.lowerBound;
	upperBound = rvalue.upperBound;
	lowerBoundMode = rvalue.lowerBoundMode;
	upperBoundMode = rvalue.upperBoundMode;
	return *this;
}

/*
 *
 */
bool Interval::contains(const Value& value) const
{
	return
		(lowerBoundMode == Unbounded or
			(lowerBound.getType() == value.getType() and	// no type checking for infinity
				((lowerBoundMode == IncludeBound and lowerBound <= value) or
				(lowerBoundMode == ExcludeBound and lowerBound < value)))) and
		(upperBoundMode == Unbounded or
			(upperBound.getType() == value.getType() and // no type checking for infinity
				((upperBoundMode == IncludeBound and value <= upperBound) or
				(upperBoundMode == ExcludeBound and value < upperBound))));
}

/*
 *
 */
bool Interval::equals(const Range& rvalue) const
{
	const Interval* interval = dynamic_cast<const Interval*>(&rvalue);
	return interval != 0 and
		lowerBoundMode == interval->lowerBoundMode and
		(lowerBoundMode == Unbounded or
			(lowerBound.getType() == interval->lowerBound.getType() and
			lowerBound == interval->lowerBound)) and
		upperBoundMode == interval->upperBoundMode and
		(upperBoundMode == Unbounded or
			(upperBound.getType() == interval->upperBound.getType() and
			upperBound == interval->upperBound));
}

/*
 *
 */
RangePtr Interval::clone() const
{
	return RangePtr(new Interval(*this));
}

/*
 *
 */
std::string Interval::toString() const
{
	std::ostringstream buffer;
	buffer << *this;
	return buffer.str();
}

/*
 *
 */
void Interval::fromString(const std::string& str)
{
	Grammar::IntervalGrammar::fromString(*this, str);
}

/*
 *
 */
void Interval::write(std::ostream& stream) const
{
	if(lowerBoundMode == Unbounded)
	{
		stream << Interval::LBOUND_EXC << MINUS_INFINITY;
	}
	else
	{
		stream
			<< (lowerBoundMode == IncludeBound ? Interval::LBOUND_INC : Interval::LBOUND_EXC)
			<< Value::DelimitedFormat(lowerBound, Interval::RESERVED_CHARACTERS);
	}
	stream << Interval::SEPARATOR << ' ';

	if(upperBoundMode == Unbounded)
	{
		stream << PLUS_INFINITY << Interval::UBOUND_EXC;
	}
	else
	{
		stream
			<< Value::DelimitedFormat(upperBound, Interval::RESERVED_CHARACTERS)
			<< (upperBoundMode == Interval::IncludeBound ? Interval::UBOUND_INC : Interval::UBOUND_EXC);
	}
}

}
