// =============================================================================
// ScalarClassInfo.cpp
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

#include "ScalarClassInfo.h"
#include "Range.h"
#include "Utils/Io/PointerValueStreamer.hpp"

using Utils::Io::writeValue;

namespace CCfg
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::ScalarClassInfo class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ScalarClassInfo::ScalarClassInfo()
	: ClassInfo(), unit(), range()
{
}

/**
 * Copy constructor.
 */
ScalarClassInfo::ScalarClassInfo(const ScalarClassInfo& rvalue)
	: ClassInfo(), unit(rvalue.unit),
	range(rvalue.range.get() != 0 ? rvalue.range->clone() : RangePtr())
{
}

/**
 * Destructor.
 */
ScalarClassInfo::~ScalarClassInfo()
{
}

/**
 * Assignment operator.
 */
ScalarClassInfo& ScalarClassInfo::operator=(const ScalarClassInfo& rvalue)
{
	unit = rvalue.unit;
	RangePtr temp = rvalue.range->clone();
	range = temp; // hack for vxworks which doesn't support assignment from return!!
	return *this;
}

/*
 *
 */
bool ScalarClassInfo::equals(const ClassInfo& rvalue) const
{
	const ScalarClassInfo* other = dynamic_cast<const ScalarClassInfo*>(&rvalue);
	return other != 0 and unit == other->unit and
		Range::equivalent(range.get(), other->range.get());
}

/*
 *
 */
std::auto_ptr<ClassInfo> ScalarClassInfo::clone() const
{
	return std::auto_ptr<ClassInfo>(new ScalarClassInfo(*this));
}

/*
 *
 */
void ScalarClassInfo::write(std::ostream& stream) const
{
	stream << "{unit=\"" << unit << "\" range=" << writeValue(range.get()) << '}';
}

}
