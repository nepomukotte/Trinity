// =============================================================================
// Scalar.h
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

#ifndef CCfg_Scalar_h_INCLUDED
#define CCfg_Scalar_h_INCLUDED

#include "Attribute.h"
#include "Value.h"
#include <string>
#include <memory>

namespace CCfg
{
///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class describes a scalar in a tree-like configuration document.
 * The scalar is the class that stores values (along with a unit) in a
 * configuration. The type of values may change, but we can make assumption
 * on the type when we retrieve the value. If the assumption proves false,
 * an exception is thrown. Alternatively, we can force conversion of the value
 * to the desired type.
 *
 * @author fchateau
 */
class Scalar : public Attribute
{
	Value val;	///< the variant value of this scalar attribute.

public:
	// @{ construction/destruction
	explicit Scalar(const std::string& name = "unnamed", const std::string& index = "");
	Scalar(const std::string& name, const std::string& index, const Value& value);
	Scalar(const Scalar& value);
	~Scalar();
	Scalar& operator =(const Scalar& value);	// assignment operator
	Scalar& operator =(const Attribute& value);	// overriden base class assignment
	std::auto_ptr<Scalar> clone(Recursivity recursivity = DEEP) const;
	// @}

	// @{ comparison
	bool equals(const Attribute& attr, MatchLevel level = IDENTICAL, Recursivity recursivity = DEEP) const;
	// @}

	// @{ visitation
	bool accept(Visitor& visitor) const;
	// @}

	// @{ getters
	const std::string& getUnit() const;
	Range* getRange() const;

	const Value&	getValue() const;
	ValueTypeId		getValueTypeId() const;
	const String&	getStringValue() const;
	Int		getIntValue() const;
	Hex		getHexValue() const;
	Real	getRealValue() const;
	Bool	getBoolValue() const;
	// @}

	// @{ conversion getters
	const String&	getValueAsString() const;
	Int		getValueAsInt() const;
	Hex		getValueAsHex() const;
	Real	getValueAsReal() const;
	Bool	getValueAsBool() const;
	// @}

	// @{ setters
	void setUnit(const std::string& value);
	void setRange(RangePtr range);
	void setRange(const Range& range);

	void setValue(const String& value);
	void setValue(const Value& value);
	void setValue(ValueTypeId type, const String& value);
	void setStringValue(const String& value);

	void setIntValue(Int value);
	void setHexValue(Hex value);
	void setRealValue(Real value);
	void setBoolValue(Bool value);

	void setIntValue(const String& value);
	void setHexValue(const String& value);
	void setRealValue(const String& value);
	void setBoolValue(const String& value);
	// @}

	// @{ assignment operators
	Scalar& operator =(const String& value);
	Scalar& operator =(Int value);
	Scalar& operator =(Hex value);
	Scalar& operator =(Real value);
	Scalar& operator =(Bool value);
	Scalar& operator =(const char* value);
	Scalar& operator =(const Value& value);
	Scalar& operator =(short value);
	Scalar& operator =(long value);
	Scalar& operator =(unsigned short value);
	Scalar& operator =(unsigned long value);
	// @}

	// @{ type information
	bool isObject() const;
	AttributeTypeId getType() const;
	// @}

private:
	Scalar* cloneImpl(Recursivity recursivity = DEEP) const;
	void updateClassInfo();
};

}

#endif // CCfg_Scalar_h_INCLUDED
