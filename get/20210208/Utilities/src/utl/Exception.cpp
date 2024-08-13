// =============================================================================
// Exception.cpp
// -----------------------------------------------------------------------------
// class utl::Exception
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, shebli.anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "Exception.h"
#include <sstream>

namespace utl
{
// =============================================================================
// utl::Exception class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
Exception::Exception()
{
}

/**
 * Copy constructor.
 * @param value the exception to copy
 */
Exception::Exception(const Exception& value)
	: std::exception(), message(value.message)
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
		std::ostringstream stream;
		writeMessage(stream);
		message = stream.str();
	}
	return message;
}

/**
 * @fn void Exception::writeMessage(std::ostream& stream) const
 * Writes an error message for the exception on a stream.
 * This method must be implemented by each subclass.
 * @param stream the stream where to write the error message.
 */

} // namespace utl


