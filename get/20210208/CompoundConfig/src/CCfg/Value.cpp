// =============================================================================
// Value.cpp
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

#include "Value.h"
#include "Grammar/ValueGrammar.h"
#include "TypeInfo.h"
#include "Utils/Math/float.h"
#include <sstream>
#include <iomanip>
#include <limits>
#include <cstring>
#include <cassert>
#include <cmath>

using std::istringstream;
using std::ostringstream;
using std::string;
using Utils::Math::equals;

namespace CCfg
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::Value class methods body
///////////////////////////////////////////////////////////////////////////////

const Real Value::EPSILON = pow(10.0, -std::numeric_limits<Real>::digits10)*0.5;
const Real Value::POS_EXP_LIMIT = pow(10.0, std::numeric_limits<Real>::digits10);
const char Value::STRING_DELIMITER;

/**
 * Default constructor.
 */
Value::Value()
	: valueType(STRING), realValue(0.0), textValue()
{

}

/**
 *  Full constructor for C-strings.
 * @param value the string value that the class will hold.
 */
Value::Value(const char* value)
	: valueType(STRING), realValue(0.0), textValue(value)
{
}

/**
 * Full constructor for strings.
 * @param value the string value that the class will hold.
 */
Value::Value(const String& value)
	: valueType(STRING), realValue(0.0), textValue(value)
{
}

/**
 * Full constructor for integers.
 * @param value the integer value that the class will hold.
 */
Value::Value(Int value)
	: valueType(STRING), realValue(0.0), textValue()
{
	setInt(value);
}

/**
 * Full constructor for hexadecimal numbers.
 * @param value the hexadecimal value that the class will hold.
 */
Value::Value(Hex value)
	: valueType(STRING), realValue(0.0), textValue()
{
	setHex(value);
}

/**
 * Full constructor for real numbers.
 * @param value the real value that the class will hold.
 */
Value::Value(Real value)
	: valueType(STRING), realValue(0.0), textValue()
{
	setReal(value);
}

/**
 * Full constructor for booleans.
 * @param value the boolean value that the class will hold.
 */
Value::Value(Bool value)
	: valueType(STRING), realValue(0.0), textValue()
{
	setBool(value);
}

/**
 * Copy constructor.
 * @param value the config Value to copy.
 */
Value::Value(const Value& value)
	: valueType(STRING), realValue(0.0), textValue()
{
	setValue(value);
}

/**
 * Destructor.
 */
Value::~Value()
{
}

/**
 * Assignement operator: overwrites this value with a copy of the specified one.
 * @param value the config Value to copy.
 */
Value& Value::operator=(const Value& value)
{
	setValue(value);
	return *this;
}

/**
 * Equality operator.
 * @return Returns true if internal values are equals, false otherwise.
 */
bool Value::operator==(const Value& rvalue) const
{
	if(valueType != rvalue.getType())
	{
		throw IncompatibleTypes(*this, rvalue);
	}
	switch(valueType)
	{
	case INTEGER:
		return intValue == rvalue.intValue;
	case HEXADECIMAL:
		return hexValue == rvalue.hexValue;
	case REAL:
		return realValue == rvalue.realValue;
	case BOOLEAN:
		return boolValue == rvalue.boolValue;
	case STRING:
		return textValue == rvalue.textValue;
	default:
		assert(false);
		return false;
	}
}

/**
 *
 */
bool Value::operator<(const Value& rvalue) const
{
	if(valueType != rvalue.getType())
	{
		throw IncompatibleTypes(*this, rvalue);
	}
	switch(valueType)
	{
	case INTEGER:
		return intValue < rvalue.intValue;
	case HEXADECIMAL:
		return hexValue < rvalue.hexValue;
	case REAL:
		return realValue < rvalue.realValue;
	case BOOLEAN:
		return boolValue < rvalue.boolValue;
	case STRING:
		return textValue < rvalue.textValue;
	default:
		assert(false);
		return false;
	}
}

/**
 * Getter of the property "valueType".
 * @return Returns the type of the value in the value (integer, real, etc...)
 */
ValueTypeId Value::getType() const
{
	return valueType;
}

/**
 * Getter of the internal value, assuming its exact type.
 * @return Returns the internal value.
 * @throws Value::UnexpectedType if the internal value is not a string.
 */
const String& Value::getString() const
{
	if(valueType != STRING)
	{
		throw UnexpectedType(*this, STRING);
	}
	return getText();
}

/**
 * Getter of the internal value, assuming its exact type.
 * @return Returns the internal value.
 * @throws Value::UnexpectedType if the internal value is not an integer.
 */
Int Value::getInt() const
{
	if(valueType != INTEGER)
	{
		throw UnexpectedType(*this, INTEGER);
	}
	return intValue;
}

/**
 * Getter of the internal value, assuming its exact type.
 * @return Returns the internal value.
 * @throws Value::UnexpectedType if the internal value is not an hexadecimal number.
 */
Hex Value::getHex() const
{
	if(valueType != HEXADECIMAL)
	{
		throw UnexpectedType(*this, HEXADECIMAL);
	}
	return hexValue;
}

/**
 * Getter of the internal value, assuming its exact type.
 * @return Returns the internal value.
 * @throws Value::UnexpectedType if the internal value is not a real number.
 */
Real Value::getReal() const
{
	if(valueType != REAL)
	{
		throw UnexpectedType(*this, REAL);
	}
	return realValue;
}

/**
 * Getter of the internal value, assuming its exact type.
 * @return Returns the internal value.
 * @throws Value::UnexpectedType if the internal value is not a boolean.
 */
bool Value::getBool() const
{
	if(valueType != BOOLEAN)
	{
		throw UnexpectedType(*this, BOOLEAN);
	}
	return boolValue;
}

/**
 * Getter the internal value, with conversion if necessary.
 * @return Returns the internal value converted to string.
 */
const String& Value::getAsString() const
{
	return getText();
}

/**
 * Getter the internal value, with conversion if necessary.
 * @return Returns the internal value converted to integer.
 */
Int Value::getAsInt() const
{
	Int intVal;
	switch(valueType)
	{
	case INTEGER:intVal = intValue;
		break;
	case HEXADECIMAL:intVal = static_cast<Int>(hexValue);
		break;
	case REAL:intVal = static_cast<Int>(realValue);
		break;
	case BOOLEAN:intVal = static_cast<Int>(boolValue);
		break;
	case STRING:
		{
			intVal = 0;
			istringstream strVal(getText());
			strVal >> intVal;
			// if conversion fails or is uncomplete
			if(strVal.fail() or not strVal.eof())
			{
				throw BadConversion(*this, INTEGER);
			}
		}
		break;
	case VALUETYPE_INVALID:
	default:
		throw BadConversion(*this, INTEGER);
	}
	return intVal;
}

/**
 * Getter the internal value, with conversion if necessary.
 * @return Returns the internal value converted to hexadecimal.
 */
Hex Value::getAsHex() const
{
	Hex hexVal;
	switch(valueType)
	{
	case INTEGER:hexVal = static_cast<Hex>(intValue);
		break;
	case HEXADECIMAL:hexVal = hexValue;
		break;
	case REAL:hexVal = static_cast<Hex>(realValue);
		break;
	case BOOLEAN:hexVal = static_cast<Hex>(boolValue);
		break;
	case STRING:
		{
			hexVal = 0u;
			istringstream strVal(getText());
			strVal >> std::hex >> hexVal;
			if(strVal.fail() or not strVal.eof())
			{
				throw BadConversion(*this, HEXADECIMAL);
			}
		}
		break;
	case VALUETYPE_INVALID:
	default:
		throw BadConversion(*this, HEXADECIMAL);
	}
	return hexVal;
}

/**
 * Getter the internal value, with conversion if necessary.
 * @return Returns the internal value converted to real.
 */
Real Value::getAsReal() const
{
	Real realVal;
	switch(valueType)
	{
	case INTEGER:realVal = static_cast<Real>(intValue);
		break;
	case HEXADECIMAL:realVal = static_cast<Real>(hexValue);
		break;
	case REAL:realVal = realValue;
		break;
	case BOOLEAN:realVal = boolValue ? 1.0 : 0.0;
		break;
	case STRING:
		{
			realVal = 0.0;
			istringstream strVal(getText());
			strVal >> realVal;
			if(strVal.fail() or not strVal.eof())
			{
				throw BadConversion(*this, REAL);
			}
		}
		break;
	case VALUETYPE_INVALID:
	default:
		throw BadConversion(*this, REAL);
	}
	return realVal;
}

/**
 * Getter the internal value, with conversion if necessary.
 * @return Returns the internal value converted to boolean.
 */
Bool Value::getAsBool() const
{
	Bool boolVal;
	switch(valueType)
	{
	case INTEGER:boolVal = (intValue != 0);
		break;
	case HEXADECIMAL:boolVal = (hexValue != 0);
		break;
	case REAL:boolVal = (realValue != 0.0);
		break;
	case BOOLEAN:boolVal = boolValue;
		break;
	case STRING:
		{
			boolVal = false;
			char leftOverChar;
			istringstream strVal(getText());
			if (not (strVal >> std::boolalpha >> boolVal) or strVal.get(leftOverChar))
			{
				throw BadConversion(*this, BOOLEAN);
			}
		}
		break;
	case VALUETYPE_INVALID:
	default:
		throw BadConversion(*this, BOOLEAN);
	}
	return boolVal;
}


/**
 * Setter of the internal value. Automatically determine the appropriate internal type.
 * @param value the new value for this property.
 */
void Value::setValue(const String& value)
{
	Grammar::ValueGrammar::fromString(*this, value);
}

/**
 * Defines the internal value with a string, while forcing its interpretation.
 */
void Value::setValue(ValueTypeId type, const String& value)
{
	switch(type)
	{
	case INTEGER:setInt(value);
		break;
	case HEXADECIMAL:setHex(value);
		break;
	case REAL:setReal(value);
		break;
	case BOOLEAN:setBool(value);
		break;
	case STRING:setString(value);
		break;
	case VALUETYPE_INVALID:
		assert(false);
		break;
	default:
		assert(false);
		break;
	}
}

/**
 * Defines the internal value and set its type to string.
 */
void Value::setString(const String& value)
{
	setText(value);
	valueType = STRING;
}

/**
 * Defines the internal value and set its type to integer.
 */
void Value::setInt(int value)
{
	intValue = value;
	valueType = INTEGER;
	// updates the string representation of the value
	ostringstream strVal;
	strVal << value;
	setText(strVal.str());
}

/**
 * Defines the internal value and set its type to hexadecimal.
 */
void Value::setHex(Hex value)
{
	hexValue = value;
	valueType = HEXADECIMAL;
	// updates the string representation of the value
	ostringstream strVal;
	strVal << "0x" << std::hex << std::uppercase << value;
	setText(strVal.str());
}

/**
 * Defines the internal value and set its type to real.
 */
void Value::setReal(Real value)
{
	realValue = value;
	valueType = REAL;
	// updates the string representation of the value
	ostringstream strVal;
	if(equals(value, floor(value), EPSILON) and fabs(value) < POS_EXP_LIMIT)
	{
		strVal << std::setprecision(1) << std::fixed << std::showpoint << value;
	}
	else
	{
		strVal << std::setprecision(std::numeric_limits<Real>::digits10) << value;
	}
	setText(strVal.str());
}

/**
 * Defines the internal value and set its type to boolean.
 */
void Value::setBool(Bool value)
{
	boolValue = value;
	valueType = BOOLEAN;
	// updates the string representation of the value
	ostringstream strVal;
	strVal << std::boolalpha << value;
	setText(strVal.str());
}

/**
 * Defines the internal value using a string interpreted as an integer.
 */
void Value::setInt(const String& value)
{
	Value stringInt(value);
	setInt(stringInt.getAsInt());
}

/**
 * Defines the internal value using a string interpreted as an hexadecimal.
 */
void Value::setHex(const String& value)
{
	Value stringHex(value);
	setHex(stringHex.getAsHex());
}

/**
 * Defines the internal value using a string interpreted as a real.
 */
void Value::setReal(const String& value)
{
	Value stringReal(value);
	setReal(stringReal.getAsReal());
}

/**
 * Defines the internal value using a string interpreted as a boolean.
 */
void Value::setBool(const String& value)
{
	Value stringBool(value);
	setBool(stringBool.getAsBool());
}

/**
 * Returns the string representation of the value.
 */
const std::string& Value::getText() const
{
	return textValue;
}

/**
 * Defines the text representation of the value.
 * @param value the new text representation.
 */
void Value::setText(const std::string& value)
{
	textValue = value;
}

/**
 * Defines the internal value of this attribute, assuming it is a value attribute.
 * @return Returns a reference to the value that has been altered.
 */
Value& Value::operator =(const String& value)
{
	setString(value);
	return *this;
}

/**
 * Defines the internal value of this attribute, assuming it is a value attribute.
 * @return Returns a reference to the value that has been altered.
 */
Value& Value::operator =(Int value)
{
	setInt(value);
	return *this;
}

/**
 * Defines the internal value of this attribute, assuming it is a value attribute.
 * @return Returns a reference to the value that has been altered.
 */
Value& Value::operator =(Hex value)
{
	setHex(value);
	return *this;
}

/**
 * Define the internal value of this attribute, assuming it is a value attribute.
 * @return Returns a reference to the value that has been altered.
 */
Value& Value::operator =(Real value)
{
	setReal(value);
	return *this;
}

/**
 * Define the internal value of this attribute, assuming it is a value attribute.
 * @return Returns a reference to the value that has been altered.
 */
Value& Value::operator =(Bool value)
{
	setBool(value);
	return *this;
}

/**
 * Defines the internal value of this attribute, assuming it is a value attribute.
 * @return Returns a reference to the value that has been altered.
 */
Value& Value::operator =(const char* value)
{
	setString(value);
	return *this;
}

/**
 * Copies a Value object.
 * It is a convenience to avoid handling different possible types each time.
 * @param value the value containing the value to copy.
 */
void Value::setValue(const Value& value)
{
	valueType = value.getType();
	switch(value.getType())
	{
	case INTEGER:intValue = value.intValue;
		break;
	case HEXADECIMAL:hexValue = value.hexValue;
		break;
	case REAL:realValue = value.realValue;
		break;
	case BOOLEAN:boolValue = value.boolValue;
		break;
	case STRING:// the next instruction: setTextValue does the job.
		break;
	case VALUETYPE_INVALID:
		assert(false);
		break;
	default:
		assert(false);
		break;
	}
	// in any case we have to copy text representation.
	setText(value.getText());
}






////////////////////////////////////////////////////////////////////////////////
// CCfg::Value::DelimitedFormat class methods body
////////////////////////////////////////////////////////////////////////////////

const char Value::DelimitedFormat::SPECIAL_CHARS[] = "\b\t\n\f\r\\\"\'";
const char* const Value::DelimitedFormat::ESCAPED_SPECIAL_CHARS[sizeof(SPECIAL_CHARS)] =
{
		"\\b",
		"\\t",
		"\\n",
		"\\f",
		"\\r",
		"\\\\",
		"\\\"",
		"\\\'"
};

/**
 * Returns whether the value to format needs delimiters.
 * @note Only string values can require delimiters.
 */
bool Value::DelimitedFormat::needDelimiters() const
{
	if(value.getType() == STRING)
	{
		const std::string& str = value.getString();
		for(size_t i = 0; i < str.length(); ++i)
		{
			if(strchr(restrictedChars, str[i]) or
				strchr(Value::DelimitedFormat::SPECIAL_CHARS, str[i]) != 0 or
				iscntrl(str[i]))
			{
				return true;
			}
		}
	}
	return false;
}

/**
 * Stream insertion operator.
 */
std::ostream& operator<<(std::ostream& stream, const Value::DelimitedFormat& format)
{
	if(format.needDelimiters())
	{
		const std::string& str = format.value.getString();
		stream << Value::STRING_DELIMITER;
		for(size_t i = 0; i < str.length(); ++i)
		{
			char c = str[i];
			const char* special = strchr(Value::DelimitedFormat::SPECIAL_CHARS, c);
			if(special != 0)
			{
				stream << Value::DelimitedFormat::ESCAPED_SPECIAL_CHARS[special - Value::DelimitedFormat::SPECIAL_CHARS];
			}
			else
			{
				if(iscntrl(c))
				{
					stream << "\\x" << std::hex << static_cast<int>(c);
				}
				else
				{
					stream << c;
				}
			}
		}
		return stream << Value::STRING_DELIMITER;
	}
	return stream << format.value;
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Value::UnexpectedType class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param value the value in fault.
 * @param expectedType the internal value type expected by the client in its call.
 */
Value::UnexpectedType::UnexpectedType(Value value, ValueTypeId expectedType)
	: value(value), expectedType(expectedType)
{
}

/**
 * Copy constructor.
 * @param value the exception to copy
 */
Value::UnexpectedType::UnexpectedType(const UnexpectedType& value)
	: Exception(value),
	value(value.value), expectedType(value.expectedType)
{
}

/**
 * Destructor.
 */
Value::UnexpectedType::~UnexpectedType() throw()
{
}

/**
 * Assignment operator.
 * @param e the exception to copy
 */
Value::UnexpectedType& Value::UnexpectedType::operator =(const UnexpectedType& e)
{
	Exception::operator =(e);
	expectedType = e.expectedType;
	value = e.value;
	return *this;
}

/**
 * Getter of the property "value".
 * @return Returns the value where type assumption proved false.
 */
const Value& Value::UnexpectedType::getValue() const
{
	return value;
}

/**
 * Getter of the property "expectedType".
 * @return Returns the expected type.
 */
ValueTypeId Value::UnexpectedType::getExpectedType() const
{
	return expectedType;
}

std::string Value::UnexpectedType::formatMessage() const
{
	ostringstream message;
	message << "Unexpected type: the value \"" << value.getAsString() << "\""
		" is of type " << TypeInfo::getDisplayName(value.getType()) <<
		" but " << TypeInfo::getDisplayName(expectedType) << " was assumed.";
	return message.str();
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Value::IncompatibleTypes class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param lvalue the left value of the operation.
 * @param rvalue the right value of the operation.
 * @param operation the type of the failed operation.
 */
Value::IncompatibleTypes::IncompatibleTypes(const Value& lvalue, const Value& rvalue,
	Operation operation)
	: operation(operation), lvalue(lvalue), rvalue(rvalue)
{
}

/**
 * Copy constructor.
 */
Value::IncompatibleTypes::IncompatibleTypes(const IncompatibleTypes& value)
	: Exception(value), lvalue(value.lvalue), rvalue(value.rvalue)
{
}

/**
 * Destructor.
 */
Value::IncompatibleTypes::~IncompatibleTypes() throw()
{
}

/**
 * Assignment operator.
 */
Value::IncompatibleTypes& Value::IncompatibleTypes::operator =(const IncompatibleTypes& value)
{
	Exception::operator=(value);
	lvalue = value.lvalue;
	rvalue = value.rvalue;
	return *this;
}

/**
 * Returns the left value of the operation.
 */
const Value& Value::IncompatibleTypes::getLValue() const
{
	return lvalue;
}

/**
 * Returns the right value of the operation.
 */
const Value& Value::IncompatibleTypes::getRValue() const
{
	return rvalue;
}

/**
 * Returns the type of the failed operation.
 */
Value::IncompatibleTypes::Operation Value::IncompatibleTypes::getOperation() const
{
	return operation;
}

/*
 *
 */
std::string Value::IncompatibleTypes::formatMessage() const
{
	std::ostringstream message;
	message << "Cannot ";
	switch(operation)
	{
	case Comparison:
		message << "compare";
		break;
	default:
		assert(false);
	}
	message << ' ' << lvalue << " with " << rvalue
		<< " because they do not have the same type.";
	return message.str();
}



///////////////////////////////////////////////////////////////////////////////
// CCfg::Value::BadConversion class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor: initializes informations about the exception.
 * @param value the value in fault.
 * @param targetType the target type in conversion.
 */
Value::BadConversion::BadConversion(Value value, ValueTypeId targetType)
	: value(value), targetType(targetType)
{

}

/**
 * Copy constructor.
 */
Value::BadConversion::BadConversion(const BadConversion& rvalue)
	: Exception(rvalue),
	value(rvalue.value), targetType(rvalue.targetType)
{
}

/**
 * Destructor.
 */
Value::BadConversion::~BadConversion() throw()
{
}

/**
 * Assignment operator.
 */
Value::BadConversion& Value::BadConversion::operator =(const Value::BadConversion& rvalue)
{
	Exception::operator =(rvalue);
	targetType = rvalue.targetType;
	value = rvalue.value;
	return *this;
}

/**
 * Getter of the property "value".
 * @return Returns the value which cannot be converted.
 */
const Value& Value::BadConversion::getValue() const
{
	return value;
}

/**
 * Getter of the property "targetType".
 * @return Returns the target type in conversion.
 */
ValueTypeId Value::BadConversion::getTargetType() const
{
	return targetType;
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string Value::BadConversion::formatMessage() const
{
	ostringstream message;
	message << "Bad conversion: "
		"the " << TypeInfo::getDisplayName(value.getType()) << " "
		"value \"" << value.getAsString() << "\" "
		"cannot be converted to " << TypeInfo::getDisplayName(targetType) << '.';
	return message.str();
}


}
