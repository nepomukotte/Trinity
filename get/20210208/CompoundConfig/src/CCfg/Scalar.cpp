// =============================================================================
// Scalar.cpp
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

#include "View/ClassIterator.h"
#include "Scalar.h"
#include "Attribute.h"
#include "Visitor.h"
#include "UnsupportedOperation.h"
#include "ScalarClassInfo.h"
#include "Range.h"

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Scalar class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor: creates a scalar with the specified name and index but no value.
 * @param name the name of the new config scalar.
 * @param index the index of the new config scalar.
 */
Scalar::Scalar(const std::string& name, const std::string& index)
	: Attribute(name, index), val()
{
}

/**
 * Full constructor.
 * @param name the name of the new config scalar.
 * @param index the index of the new config scalar.
 * @param value the value of the new config scalar.
 */
Scalar::Scalar(const std::string& name, const std::string& index, const Value& value)
	: Attribute(name, index), val(value)
{
}

/**
 * Copy constructor: makes a copy of a scalar.
 * @param scalar the config scalar to copy.
 */
Scalar::Scalar(const Scalar& scalar)
	: Attribute(scalar), val(scalar.val)
{
}

/**
 * Destructor.
 */
Scalar::~Scalar()
{
}

/**
 * Assignement operator.
 * @param scalar the config scalar to copy.
 */
Scalar& Scalar::operator =(const Scalar& scalar)
{
	Attribute::operator =(scalar);
	val = scalar.val;
	return *this;
}

/**
 * Assignement operator from a scalar treated as an attribute.
 * @param attr the config scalar to copy.
 * @throws UnsupportedOperation if the right value is not a CCfg::Scalar.
 */
Scalar& Scalar::operator =(const Attribute& attr)
{
	Attribute::operator =(attr);
	val = attr.getValue();
	return *this;
}

/**
 * Makes a clone of this scalar.
 * @note
 * The clone is not included as a child of the parent of its model.
 * @param recursivity it has no effect on a scalar, because recursivity ends here.
 */
Scalar* Scalar::cloneImpl(Recursivity /*recursivity*/) const
{
	return new Scalar(*this);
}

/*
 *
 */
void Scalar::updateClassInfo()
{
	View::ClassIterator it(*this);
	Attribute* otherClassAttr = 0;
	if(it.hasNext()) // find the first scalar belonging to the same "equivalence class".
	{
		it.next();
		otherClassAttr = it.get(); // fetch it
		if(otherClassAttr == this) // if the scalar found itself
		{
			if(it.hasNext()) // try to find another one
			{
				it.next();
				otherClassAttr = it.get(); // fetch it
			}
			else
			{
				otherClassAttr = 0; // cannot find one
			}
		}
	}

	if(otherClassAttr != 0) // if we found another scalar belonging to the same class
	{
		// retrieve its class info (that should be shared by all instances of this class).
		// it can be null.
		const ClassInfoPtr& commonClassInfo = otherClassAttr->getClassInfo();

		if(classInfo != commonClassInfo) // if our class info is different from that one.
		{
			// if we don't have a class info or if both are not null, and ours
			// has the same content as the one we found.
			if(classInfo.get() == 0 or
				(commonClassInfo.get() != 0 and classInfo->equals(*commonClassInfo)))
			{
				// set our class info so that it becomes shared with other class members.
				classInfo = commonClassInfo;
			}
			else
			{
				// otherwise we have a class info that is different from the already
				// existing class info for this class.
				throw IncompatibleClassInfo(this, otherClassAttr);
			}
		} // else: nothing to do, the class info is fine.
	}
	else // this is the only member of its class.
	{
		if(classInfo.get() != 0) // if it has a class info
		{
			if(not classInfo.unique()) // if the class info is shared
			{
				// in case the class info is shared, it means that the scalar
				// still have a link to a class info that may come from another
				// configuration tree.
				// we don't want to keep a dependency with another configuration:
				// the link must be severed, so we make a copy.
				std::auto_ptr<ClassInfo> classInfoCopy = classInfo->clone();
				classInfo = classInfoCopy;
			}
		}
	}
}

/**
 * Makes a clone of this scalar.
 * @note
 * The clone is not included as a child of the parent of its model.
 * @param recursivity it has no effect on a scalar, because recursivity ends here.
 * @return Returns a newly allocated copy of this scalar.
 */
std::auto_ptr<Scalar> Scalar::clone(Recursivity recursivity) const
{
	return std::auto_ptr<Scalar>(cloneImpl(recursivity));
}

/*
 * (see inherited doc)
 */
bool Scalar::equals(const Attribute& attr, MatchLevel level, Recursivity recursivity) const
{
	const Scalar* scalar = dynamic_cast<const Scalar*>(&attr);
	if(scalar == 0) // check if attr is also a Scalar.
	{
		return false;
	}
	const ClassInfoPtr& leftClass = getClassInfo();
	const ClassInfoPtr& rightClass = attr.getClassInfo();
	return Attribute::equals(attr, level, recursivity) and
		getValueTypeId() == scalar->getValueTypeId() and // compare type before comparing values because comparing values of different types returns an exception.
		val == scalar->getValue() and
		((leftClass.get() == 0 and rightClass.get() == 0) or
		(leftClass.get() != 0 and leftClass->equals(*rightClass)));
}

/**
 * Implement visitor support.
 * Call the visitor method corresponding to scalars.
 * @param visitor the visitor object.
 * @return Returns false if the visit was interrupted by the visitor (if
 * the visitor method returned false)
 */
bool Scalar::accept(Visitor& visitor) const
{
	return visitor.onScalar(*this);
}

/**
 * Getter of the property "unit".
 * @return Returns the physical unit of the value.
 */
const std::string& Scalar::getUnit() const
{
	static const std::string defaultUnit;
	if(classInfo != 0)
	{
		return static_cast<ScalarClassInfo*>(classInfo.get())->getUnit();
	}
	else
	{
		return defaultUnit;
	}
}

/*
 * Returns the range of the scalar value.
 */
Range* Scalar::getRange() const
{
	if(classInfo != 0)
	{
		return static_cast<ScalarClassInfo*>(classInfo.get())->getRange();
	}
	else
	{
		return 0;
	}
}

/**
 * Getter of the property "val".
 * @return Returns the scalar value.
 */
const Value& Scalar::getValue() const
{
	return val;
}

/**
 * Getter of the property "valueType".
 * @return Returns the type of the value in the scalar (integer, real, etc...)
 */
ValueTypeId Scalar::getValueTypeId() const
{
	return val.getType();
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a string.
 */
const String& Scalar::getStringValue() const
{
	return val.getString();
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an integer.
 */
Int Scalar::getIntValue() const
{
	return val.getInt();
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not an hexadecimal number.
 */
Hex Scalar::getHexValue() const
{
	return val.getHex();
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a real number.
 */
Real Scalar::getRealValue() const
{
	return val.getReal();
}

/**
 * Getter of the scalar value, assuming its exact type.
 * @return Returns the scalar value.
 * @throws Scalar::UnexpectedType if the scalar value is not a boolean.
 */
bool Scalar::getBoolValue() const
{
	return val.getBool();
}

/**
 * Getter the scalar value, with conversion if necessary.
 * @return Returns the scalar value converted to string.
 */
const String& Scalar::getValueAsString() const
{
	return val.getAsString();
}

/**
 * Getter the scalar value, with conversion if necessary.
 * @return Returns the scalar value converted to integer.
 */
Int Scalar::getValueAsInt() const
{
	return val.getAsInt();
}

/**
 * Getter the scalar value, with conversion if necessary.
 * @return Returns the scalar value converted to hexadecimal.
 */
Hex Scalar::getValueAsHex() const
{
	return val.getAsHex();
}

/**
 * Getter the scalar value, with conversion if necessary.
 * @return Returns the scalar value converted to real.
 */
Real Scalar::getValueAsReal() const
{
	return val.getAsReal();
}

/**
 * Getter the scalar value, with conversion if necessary.
 * @return Returns the scalar value converted to boolean.
 */
Bool Scalar::getValueAsBool() const
{
	return val.getAsBool();
}

/**
 * Defines the unit of the scalars sharing the same class.
 */
void Scalar::setUnit(const std::string& value)
{
	if(classInfo == 0)
	{
		setClassInfo(ClassInfoPtr(new ScalarClassInfo));
	}
	static_cast<ScalarClassInfo*>(classInfo.get())->setUnit(value);
}

/**
 * Defines the range of the equivalence class of this Scalar.
 */
void Scalar::setRange(RangePtr range)
{
	if(classInfo == 0)
	{
		setClassInfo(ClassInfoPtr(new ScalarClassInfo));
	}
	static_cast<ScalarClassInfo*>(classInfo.get())->setRange(range);
}

/**
 * Defines the range of the equivalence class of this Scalar.
 */
void Scalar::setRange(const Range& range)
{
	RangePtr rangePtr(range.clone());
	setRange(rangePtr);
}

/**
 * Setter of the scalar value. Automatically determine the appropriate internal type.
 * @param value the new value for this property.
 */
void Scalar::setValue(const String& value)
{
	val.setValue(value);
}

/**
 * Setter of the scalar value.
 * @param value the new value for this property.
 */
void Scalar::setValue(const Value& value)
{
	val.setValue(value);
}

/**
 * Setter of the scalar value. Forces type interpretation.
 */
void Scalar::setValue(ValueTypeId type, const String& value)
{
	val.setValue(type, value);
}

/**
 * Defines the scalar value and set its type to string.
 * @param value the new value for the scalar.
 */
void Scalar::setStringValue(const String& value)
{
	val.setString(value);
}

/**
 * Defines the scalar value and set its type to integer.
 * @param value the new value for the scalar.
 */
void Scalar::setIntValue(int value)
{
	val.setInt(value);
}

/**
 * Defines the scalar value and set its type to hexadecimal.
 * @param value the new value for the scalar.
 */
void Scalar::setHexValue(Hex value)
{
	val.setHex(value);
}

/**
 * Defines the scalar value and set its type to real.
 * @param value the new value for the scalar.
 */
void Scalar::setRealValue(Real value)
{
	val.setReal(value);
}

/**
 * Defines the scalar value and set its type to boolean.
 * @param value the new value for the scalar.
 */
void Scalar::setBoolValue(bool value)
{
	val.setBool(value);
}

void Scalar::setIntValue(const String& value)
{
	val.setInt(value);
}

void Scalar::setHexValue(const String& value)
{
	val.setHex(value);
}

void Scalar::setRealValue(const String& value)
{
	val.setReal(value);
}

void Scalar::setBoolValue(const String& value)
{
	val.setBool(value);
}

/**
 * Defines the scalar value of this attribute, assuming it is a scalar attribute.
 * @return Returns a reference to the scalar that has been altered.
 */
Scalar& Scalar::operator =(const String& value)
{
	val = value;
	return *this;
}

/**
 * Defines the scalar value of this attribute, assuming it is a scalar attribute.
 * @return Returns a reference to the scalar that has been altered.
 */
Scalar& Scalar::operator =(Int value)
{
	val = value;
	return *this;
}

/**
 * Defines the scalar value of this attribute, assuming it is a scalar attribute.
 * @return Returns a reference to the scalar that has been altered.
 */
Scalar& Scalar::operator =(Hex value)
{
	val = value;
	return *this;
}

/**
 * Define the scalar value of this attribute, assuming it is a scalar attribute.
 * @return Returns a reference to the scalar that has been altered.
 */
Scalar& Scalar::operator =(Real value)
{
	val = value;
	return *this;
}

/**
 * Define the scalar value of this attribute, assuming it is a scalar attribute.
 * @return Returns a reference to the scalar that has been altered.
 */
Scalar& Scalar::operator =(Bool value)
{
	val = value;
	return *this;
}

/**
 * Defines the scalar value of this attribute, assuming it is a scalar attribute.
 * @return Returns a reference to the scalar that has been altered.
 */
Scalar& Scalar::operator =(const char* value)
{
	val = value;
	return *this;
}

Scalar& Scalar::operator =(const Value& value)
{
	val = value;
	return *this;
}

Scalar& Scalar::operator =(short value)
{
	val.setInt(value);
	return *this;
}

Scalar& Scalar::operator =(long value)
{
	val.setInt(value);
	return *this;
}

Scalar& Scalar::operator =(unsigned short value)
{
	val.setHex(value);
	return *this;
}

Scalar& Scalar::operator =(unsigned long value)
{
	val.setHex(value);
	return *this;
}

/*
 * (see inherited doc)
 */
bool Scalar::isObject() const
{
	return false;
}

/*
 * (see inherited doc)
 */
AttributeTypeId Scalar::getType() const
{
	return TypeInfo::getAttributeTypeId(getValueTypeId());
}

}

