// =============================================================================
// Enum.cpp
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

#include "Enum.h"
#include "Grammar/EnumGrammar.h"
#include <sstream>

namespace CCfg
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Enum class methods body
////////////////////////////////////////////////////////////////////////////////

const char Enum::OPEN_CHAR;
const char Enum::CLOSE_CHAR;
const char Enum::SEPARATOR;
const char Enum::RESERVED_CHARACTERS[] = ",{} ";

/**
 * Default constructor.
 */
Enum::Enum()
	: values()
{
}

/**
 * Copy constructor.
 */
Enum::Enum(const Enum& rvalue)
	: Range(), values(rvalue.values)
{
}


/**
 * Destructor.
 */
Enum::~Enum()
{
}

/**
 * Assignment operator.
 */
Enum& Enum::operator=(const Enum& rvalue)
{
	values = rvalue.values;
	return *this;
}

/*
 *
 */
bool Enum::contains(const Value& value) const
{
	return values.find(value) != values.end();
}

/*
 *
 */
bool Enum::equals(const Range& rvalue) const
{
	const Enum* enumeration = dynamic_cast<const Enum*>(&rvalue);
	return enumeration != 0 and values == enumeration->getValues();
}

/*
 *
 */
RangePtr Enum::clone() const
{
	return RangePtr(new Enum(*this));
}

/*
 *
 */
std::string Enum::toString() const
{
	std::ostringstream buffer;
	buffer << *this;
	return buffer.str();
}

/*
 *
 */
void Enum::fromString(const std::string& str)
{
	Grammar::EnumGrammar::fromString(*this, str);
}

/*
 *
 */
void Enum::write(std::ostream& stream) const
{
	stream << Enum::OPEN_CHAR;
	int valueCount = 0;
	for(Set::const_iterator val = values.begin(); val != values.end(); ++val, ++valueCount)
	{
		if(valueCount != 0)
		{
			stream << Enum::SEPARATOR << ' ';
		}
		stream << Value::DelimitedFormat(*val, Enum::RESERVED_CHARACTERS);
	}
	stream << Enum::CLOSE_CHAR;
}

/**
 * Checks if the container would have homogenous type if the specified value was
 * inserted and throws an exception if it would not.
 * @throws Enum::checkTypeHomogeneity
 */
void Enum::checkTypeHomogeneity(const Value& v)
{
	if(values.size() > 0)
	{
		if(values.begin()->getType() != v.getType())
		{
			throw IncompatibleType(*this, v);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Enum::IncompatibleType class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
Enum::IncompatibleType::IncompatibleType(const Enum& reporterEnum, const Value& invalidValue)
	: reporterEnum(reporterEnum), invalidValue(invalidValue)
{
}

/**
 * Copy constructor.
 */
Enum::IncompatibleType::IncompatibleType(const IncompatibleType& rvalue)
	: Exception(rvalue)
{
	reporterEnum = rvalue.reporterEnum;
	invalidValue = rvalue.invalidValue;
}


/**
 * Destructor.
 */
Enum::IncompatibleType::~IncompatibleType() throw()
{
}

/**
 * Returns a copy of the enum that reported the error.
 */
const Enum& Enum::IncompatibleType::getReporterEnum() const
{
	return reporterEnum;
}

/**
 * Returns the type of values in the enum.
 */
ValueTypeId Enum::IncompatibleType::getEnumType() const
{
	return reporterEnum.getValues().begin()->getType();
}

/**
 * Returns the value that could not be inserted in the enum because of its type.
 */
const Value& Enum::IncompatibleType::getInvalidValue() const
{
	return invalidValue;
}

/*
 *
 */
std::string Enum::IncompatibleType::formatMessage() const
{
	std::ostringstream message;
	message << "There was an attempt to insert the value: "
		<< Value::DelimitedFormat(invalidValue, " ")
		<< " ("	<< TypeInfo::getDisplayName(getInvalidValue().getType()) << ')'
		<< " in the enum " << reporterEnum
		<< " of type \"" << TypeInfo::getDisplayName(getEnumType()) << "\".";
	return message.str();
}


}
