// =============================================================================
// ScalarClassInfo.h
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

#ifndef CCfg_ScalarClassInfo_h_INCLUDED
#define CCfg_ScalarClassInfo_h_INCLUDED

#include "ClassInfo.h"
#include "Range.h"
#include <string>
#include <memory>

namespace CCfg
{
// forward declarations
class Range;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class describe metadata shared by members of a same equivalence class of
 * Scalar.
 *
 * @author fchateau
 */
class ScalarClassInfo : public ClassInfo
{
	std::string unit;	///< the unit of all scalars of the equivalence class.
	RangePtr range;	///< the range of authorized values for all scalars of the equivalence class.

public:
	// @{ construction/destruction
	ScalarClassInfo();
	ScalarClassInfo(const ScalarClassInfo& rvalue);
	virtual ~ScalarClassInfo();
	ScalarClassInfo& operator=(const ScalarClassInfo& rvalue);
	// @}

	// @{ overriden methods
	bool equals(const ClassInfo& rvalue) const;
	std::auto_ptr<ClassInfo> clone() const;
	void write(std::ostream& stream) const;
	// @}

	// @{ accessors
	const std::string& getUnit() const;
	Range* getRange() const;
	void setUnit(const std::string& newUnit);
	void setRange(RangePtr newRange);
	// @}

	RangePtr resetRange();
};

////////////////////////////////////////////////////////////////////////////////
// CCfg::ScalarClassInfo inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Returns the unit of all scalars of the equivalence class.
 */
inline const std::string& ScalarClassInfo::getUnit() const
{
	return unit;
}

/**
 * Returns the range of authorized values for all scalars of the equivalence class.
 */
inline Range* ScalarClassInfo::getRange() const
{
	return range.get();
}

/**
 * Defines the unit of all scalars of the equivalence class.
 */
inline void ScalarClassInfo::setUnit(const std::string& newUnit)
{
	unit = newUnit;
}

/**
 * Defines the range of authorized values for all scalars of the equivalence class.
 */
inline void ScalarClassInfo::setRange(RangePtr newRange)
{
	range = newRange;
}

/**
 * Resets the range and returns its previous value.
 */
inline RangePtr ScalarClassInfo::resetRange()
{
	return range;
}

}

#endif // CCfg_ScalarClassInfo_h_INCLUDED
