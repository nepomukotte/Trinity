// =============================================================================
// AttributeClass.h
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

#ifndef CCfg_AttributeClass_INCLUDED
#define CCfg_AttributeClass_INCLUDED

#include "ValueRange.h"
#include "CCfg/TypeInfo.h"
#include <string>

namespace CCfg
{
//typedef
typedef std::auto_ptr < ValueRange > ValueRangePtr;

//previous declaration
class Attribute;


class AttributeClass
{
private:
	std::string path; ///< Qualified name
	std::string name; ///< local name
	std::string unit;
	AttributeTypeId attrTypeId; ///< OBJECT or SCALAR see TypeInfo.h
	ValueRangePtr valueRange; ///< pointers on different Ranges than the class Range

public:
	AttributeClass();
	explicit AttributeClass(const Attribute* attri); ///constructor

	///@name Getters
	//@{
	std::string getPath() const { return path; }
	std::string getName() const { return name; }
	std::string getUnit() const { return unit; }
	AttributeTypeId getTypeId() const { return attrTypeId; }
	const ValueRange* getValueRange() const { return valueRange.get(); }
    ValueRange* getValueRange() { return valueRange.get(); }
	ValueRangePtr getValueRangePtr()  { return valueRange; }
	//@}

	///@name Setters
	//@{
	void setPath(std::string s);
	void setName(std::string s);
	void setUnit(std::string s);
	void setTypeId(AttributeTypeId a);
	void setValueRange(ValueRangePtr r);
	//@}
};

///////////////////////////////////////////////////////////////////////////////
// CCfg::AttributeClass inline methods body
///////////////////////////////////////////////////////////////////////////////


inline void AttributeClass::setPath(std::string s)
{
	path = s;
}

inline void AttributeClass::setName(std::string s)
{
	name = s;
}

inline void AttributeClass::setUnit(std::string s)
{
	unit = s;
}

inline void AttributeClass::setTypeId(AttributeTypeId a)
{
	attrTypeId = a;
}

inline void AttributeClass::setValueRange(ValueRangePtr r)
{
	valueRange = r;
}


}//namespace CCfg



#endif /* CCfg_AttributeClass_INCLUDED */
