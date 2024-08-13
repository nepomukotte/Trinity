// =============================================================================
// Value.h
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

#ifndef CCfg_Value_h_INCLUDED
#define CCfg_Value_h_INCLUDED

#include "Exception.h"
#include "TypeInfo.h"
#include <boost/operators.hpp>
#include <string>
#include <ios>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

class Value : boost::totally_ordered<Value>
{
	ValueTypeId valueType;	///< the type of the current value
	/** Stores non-string values, for fast access.
	 */
	union
	{
		Int intValue;
		Hex hexValue;
		Real realValue;
		Bool boolValue;
	};
	std::string textValue;	///< always store the text representation of the current value.

	static const Real EPSILON; // the minimum difference between two numbers to consider them unequal.
	static const Real POS_EXP_LIMIT; // the smallest value requiring a positive exponent in string representation.

public:
	static const char STRING_DELIMITER = '"';

	// inner classes
	struct DelimitedFormat;
	class UnexpectedType;
	class IncompatibleTypes;
	class BadConversion;

	// construction/destruction
	Value();
	explicit Value(const char* value);
	explicit Value(const String& value);
	explicit Value(Int value);
	explicit Value(Hex value);
	explicit Value(Real value);
	explicit Value(Bool value);
	Value(const Value& value);
	~Value();
	Value& operator=(const Value& value);

	// comparison operators
	bool operator==(const Value& rvalue) const;
	bool operator<(const Value& rvalue) const;
	// other comparison operators are provided by boost::totally_ordered

	// getters
	ValueTypeId getType() const;
	const String&	getString() const;
	Int		getInt() const;
	Hex		getHex() const;
	Real	getReal() const;
	Bool	getBool() const;

	// conversion getters
	const String&	getAsString() const;
	Int		getAsInt() const;
	Hex		getAsHex() const;
	Real	getAsReal() const;
	Bool	getAsBool() const;

	// conversion operators
	operator const String&() const;
	operator Int() const;
	operator Hex() const;
	operator Real() const;
	operator Bool() const;
	operator const char*() const;

	// setters
	void setValue(const String& value);
	void setValue(const Value& value);
	void setValue(ValueTypeId type, const String& value);
	void setString(const String& value);

	void setInt(Int value);
	void setHex(Hex value);
	void setReal(Real value);
	void setBool(Bool value);

	void setInt(const String& value);
	void setHex(const String& value);
	void setReal(const String& value);
	void setBool(const String& value);

	// assignment operators
	Value& operator =(const String& value);
	Value& operator =(Int value);
	Value& operator =(Hex value);
	Value& operator =(Real value);
	Value& operator =(Bool value);
	Value& operator =(const char* value);

private:
	const std::string& getText() const;
	void setText(const std::string& value);
};

std::ostream& operator<<(std::ostream& stream, const Value& value);


/**
 * This class is a stream manipulator to output a Value in its delimited form if
 * necessary.
 * The delimited form must be used when a Value is written to stream with other
 * objects, characters, and will have to be parsed later.
 * A delimited string Value will be written between quotes if it contains at least
 * one quote or one of the characters specified in the constructor.
 */
struct Value::DelimitedFormat
{
	static const char SPECIAL_CHARS[];
	static const char* const ESCAPED_SPECIAL_CHARS[];

	const Value& value;
	const char* restrictedChars;

	DelimitedFormat(const Value& value, const char* restrictedChars = "");
	bool needDelimiters() const;
};

std::ostream& operator<<(std::ostream& stream, const Value::DelimitedFormat& format);

/**
 * This exception reports a wrong assumption about the type of the internal
 * value of a CCfg::Value.
 *
 * @author fchateau
 */
class Value::UnexpectedType : public Exception
{
	Value value;	///< the value where type assumption proved false.
	ValueTypeId expectedType;	///< the expected type.

public:
	// construction/destruction
	UnexpectedType(Value value, ValueTypeId expectedType);
	UnexpectedType(const UnexpectedType& rvalue);
	virtual ~UnexpectedType() throw();
	UnexpectedType& operator =(const UnexpectedType& rvalue);

	// getters
	const Value& getValue() const;
	ValueTypeId getExpectedType() const;

private:
	std::string formatMessage() const;
};

/**
 * This exception reports that an operation between 2 values cannot be carried out
 * because their type are different.
 */
class Value::IncompatibleTypes : public Exception
{
public:
	// typedef
	enum Operation {Comparison};
	// @{ construction/destruction
	IncompatibleTypes(const Value& lvalue, const Value& rvalue, Operation operation = Comparison);
	IncompatibleTypes(const IncompatibleTypes& rvalue);
	virtual ~IncompatibleTypes() throw();
	IncompatibleTypes& operator =(const IncompatibleTypes& rvalue);
	// @}

	// @{ getters
	const Value& getLValue() const;
	const Value& getRValue() const;
	Operation getOperation() const;
	// @}
private:
	std::string formatMessage() const;

	Operation operation; ///< the type of the failed operation.
	Value lvalue;	///< the left value of the operation.
	Value rvalue;	///< the right value of the operation.
};

/**
 * This exception reports that a type conversion failed.
 * A bad conversion happens when a some data type and/or value is completely
 * incompatible with the target type (example: converting a string "Hello"
 * to a real number) or when the conversion is uncomplete (example: converting
 * a string "640x480" to int is uncomplete).
 *
 */
class Value::BadConversion : public Exception
{
	Value value;	///< the value which cannot be converted.
	ValueTypeId targetType;	///< the target type in conversion.

public:
	// construction/destruction
	BadConversion(Value value, ValueTypeId targetType);
	BadConversion(const BadConversion& rvalue);
	virtual ~BadConversion() throw();
	BadConversion& operator =(const BadConversion& rvalue);

	// getters
	const Value& getValue() const;
	ValueTypeId getTargetType() const;

private:
	std::string formatMessage() const;
};

////////////////////////////////////////////////////////////////////////////////
// CCfg::Value inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Stream insertion operator.
 */
inline std::ostream& operator<<(std::ostream& stream, const Value& value)
{
	return stream << value.getAsString();
}

/**
 * String conversion operator.
 * @see getString()
 */
inline Value::operator const String&() const
{
	return getString();
}

/**
 * Integer conversion operator.
 * @see getInt()
 */
inline Value::operator Int() const
{
	return getInt();
}

/**
 * Hexadecimal conversion operator.
 * @see getHex()
 */
inline Value::operator Hex() const
{
	return getHex();
}

/**
 * Real conversion operator.
 * @see getReal()
 */
inline Value::operator Real() const
{
	return getReal();
}

/**
 * Boolean conversion operator.
 * @see getBool()
 */
inline Value::operator Bool() const
{
	return getBool();
}

/**
 * C-String conversion operator.
 */
inline Value::operator const char*() const
{
	return getString().c_str();
}


////////////////////////////////////////////////////////////////////////////////
// CCfg::Value::DelimitedFormat inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param value the Value to write in delimited format.
 * @param restrictedChars if the Value is a string and contains at least one of
 * characters in this string, or a quote, the string output will be surrounded
 * by quotes, and quotes will be escaped.
 */
inline Value::DelimitedFormat::DelimitedFormat(const Value& value, const char* restrictedChars)
	: value(value), restrictedChars(restrictedChars)
{
}

}

#endif // CCfg_Value_h_INCLUDED
