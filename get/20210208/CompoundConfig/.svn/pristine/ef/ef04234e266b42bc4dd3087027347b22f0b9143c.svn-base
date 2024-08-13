// =============================================================================
// ValueRange.h
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Kirell Benzi
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef CCfg_ValueRange_INCLUDED
#define CCfg_ValueRange_INCLUDED

#include "CCfg/TypeInfo.h"
#include "CCfg/Interval.h"
#include "CCfg/Enum.h"
#include <set>

namespace CCfg
{
/**
 * @enum CCfg::RangeCode
 * Different range code :
 */
enum RangeCode
{
	NO_RANGE, //no range
	INTERVAL,
	ENUM
};

/**
 * Struct used for ValueRangeInterval
 * @tparam T Type of the interval (Int, Hex, etc.)
 */
template <typename T>
struct Boundary
{
	Interval::BoundMode boundMode;
	T value;
};

/**
 * Main class for a range subdivided in 2 template classes :
 * - ValueRangeInterval<T>
 * - ValueRangeEnum<T>
 */
class ValueRange
{
public:
	ValueRange(RangeCode code, ValueTypeId valueType) : code(code), valueType(valueType) {}
	virtual ~ValueRange() {}

public:
	const RangeCode code;
	const ValueTypeId valueType; ///INTEGER, HEXADECIMAL, REAL, BOOLEAN, STRING, VALUETYPE_INVALID
};


/**
 * Interval of the AttributeClass.
 * @tparam T type for the interval
 */
template <typename T>
class ValueRangeInterval : public ValueRange
{
public:
	ValueRangeInterval(Interval* interv);
	ValueRangeInterval();
	Boundary<T>& operator [] (int index) { return bounds[index]; }
	const Boundary<T>& operator [] (int index) const { return bounds[index]; }

public:
	Boundary<T> bounds[2];
};

/**
 * Enum of the AttributeClass.
 * @tparam T type for the enum
 */
template <typename T>
class ValueRangeEnum : public ValueRange
{
public:
	ValueRangeEnum(Enum* enumRange);
	ValueRangeEnum();
	T operator [] (int index) const { return values[index]; }

public:
	std::set<T> values;
};

//String specialization
template<>
inline ValueRangeInterval<String>::ValueRangeInterval()
	: ValueRange(INTERVAL, Id<String>::val)
{
}

//Empty constructor
template <typename T>
ValueRangeInterval<T>::ValueRangeInterval()
	: ValueRange(INTERVAL, Id<T>::val)
{

}
//String specialization
template<>
inline ValueRangeInterval<String>::ValueRangeInterval(Interval* interv)
	: ValueRange(INTERVAL, Id<String>::val)
{
	bounds[0].boundMode = interv->getLowerBoundMode();
	bounds[1].boundMode = interv->getUpperBoundMode();
	bounds[0].value = interv->getLowerBound().getString();
	bounds[1].value = interv->getUpperBound().getString();
}

/**
 * Constructor
 * @param interv where the ValueRangeInterval is built from
 */
template <typename T>
ValueRangeInterval<T>::ValueRangeInterval(Interval* interv)
	: ValueRange(INTERVAL, Id<T>::val)
{
	bounds[0].boundMode = interv->getLowerBoundMode();
	bounds[1].boundMode = interv->getUpperBoundMode();

	if(interv->getLowerBoundMode() != Interval::Unbounded)
	{
		bounds[0].value = static_cast<T>(interv->getLowerBound());
	}
	if(interv->getUpperBoundMode() != Interval::Unbounded)
	{
		bounds[1].value = static_cast<T>(interv->getUpperBound());
	}

}


//String specialization
template <>
inline ValueRangeEnum<String>::ValueRangeEnum()
	: ValueRange(ENUM, Id<String>::val)
{
}

//Empty constructor
template <typename T>
ValueRangeEnum<T>::ValueRangeEnum()
	: ValueRange(ENUM, Id<T>::val)
{
}

//String specialization
template <>
inline ValueRangeEnum<String>::ValueRangeEnum(Enum* enumRange)
	: ValueRange(ENUM, Id<String>::val)
{
	const Enum::Set& valuesSet = enumRange->getValues();
	for (Enum::Set::const_iterator it = valuesSet.begin(); it != valuesSet.end(); ++it)
	{
		values.insert(it->getString());
	}
}

/**
 * Constructor
 * @param enumRange  where the ValueRangeEnum is built from
 */
template <typename T>
ValueRangeEnum<T>::ValueRangeEnum(Enum* enumRange)
	: ValueRange(ENUM, Id<T>::val)
{
	const Enum::Set& valuesSet = enumRange->getValues();
	for (Enum::Set::const_iterator it = valuesSet.begin(); it != valuesSet.end(); ++it)
	{
		values.insert(static_cast<T>(*it));
	}
}

} // namespace CCfg

#endif /* CCfg_ValueRange_INCLUDED */
