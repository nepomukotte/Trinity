// =============================================================================
// ClassInfo.h
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

#ifndef CCfg_ClassInfo_h_INCLUDED
#define CCfg_ClassInfo_h_INCLUDED

#include <boost/operators.hpp>
#include <memory>

namespace CCfg
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is the base class for metadata classes shared by members of equivalence
 * classes of Attributes.
 * These equivalence classes are based on the generic path of the attributes:
 * attributes that have the same generic path and belong to the same
 * configuration are members of the same equivalence class.
 *
 * @author fchateau
 */
class ClassInfo : boost::equality_comparable<ClassInfo>
{

public:
	ClassInfo();
	virtual ~ClassInfo();

	virtual bool equals(const ClassInfo& rvalue) const = 0;
	virtual std::auto_ptr<ClassInfo> clone() const = 0;
	virtual void write(std::ostream& stream) const = 0;

	bool operator==(const ClassInfo& rvalue) const;
};

std::ostream& operator<<(std::ostream& stream, const ClassInfo& classInfo);

////////////////////////////////////////////////////////////////////////////////
// CCfg::ClassInfo inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Equality comparison operator.
 */
inline bool ClassInfo::operator==(const ClassInfo& rvalue) const
{
	return equals(rvalue);
}

/**
 * Stream insertion operator.
 */
inline std::ostream& operator<<(std::ostream& stream, const ClassInfo& classInfo)
{
	classInfo.write(stream);
	return stream;
}

}

#endif // CCfg_ClassInfo_h_INCLUDED
