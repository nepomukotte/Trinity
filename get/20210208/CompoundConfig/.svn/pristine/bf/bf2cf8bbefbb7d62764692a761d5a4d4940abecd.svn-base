// =============================================================================
// Exception.cpp
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

#include "Exception.h"

namespace CCfg
{

/**
 * Default constructor: initializes an exception with no message.
 */
Exception::Exception()
{
}

/**
 * Implicit constructor: initializes an exception with the specified message.
 * @param message the initial message.
 */
Exception::Exception(const std::string& message) :
	message(message)
{
}

/**
 * Copy constructor.
 * @param value the exception to copy
 */
Exception::Exception(const Exception& value)
	: std::exception(value), message(value.message)
{
}

/**
 * Assignment operator.
 * @param value the exception to copy
 */
Exception& Exception::operator =(const Exception& value)
{
	message = value.message;
	return *this;
}

/**
 * Destructor.
 */
Exception::~Exception() throw()
{
}

/**
 * Returns an error message describing the exception.
 * @return Returns the error message.
 */
const char* Exception::what() const throw()
{
	return getMessage().c_str();
}

/**
 * Getter of the property "message".
 * Lazily initializes the message attribute with function formatMessage.
 * @return Returns an error message describing the cause of the exception.
 */
const std::string& Exception::getMessage() const
{
	if(message.empty())
	{
		message = formatMessage();
	}
	return message;
}

}

