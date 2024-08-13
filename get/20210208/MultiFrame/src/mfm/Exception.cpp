/* =============================================================================
 * $Id: Exception.cpp 747 2012-06-19 15:00:42Z psizun $
 * @file Exception.cpp
 * @author sizun
 * @date 15 déc. 2011
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#include "Exception.h"
#include <sstream>

namespace mfm {
//______________________________________________________________________
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
Exception::Exception(const Exception& value) :
	std::exception(value), message(value.message)
{
}

/**
 * Assignment operator.
 * @param value The exception to copy.
 * @return Returns a reference to this.
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

std::string Exception::formatMessage() const
{
	return message;
}
//______________________________________________________________________
/**
 * Constructor.
 */
ItemNotFound::ItemNotFound(size_t indexNotFound, size_t itemCount)
	: indexNotFound_(indexNotFound), itemCount_(itemCount)
{
	;
}

/**
 * Destructor.
 */
ItemNotFound::~ItemNotFound() throw()
{
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string ItemNotFound::formatMessage() const
{
	std::ostringstream message;
	message << "Item # " << indexNotFound_ << " not found. The frame only contains " << itemCount_ << "items.";
	return message.str();
}
//______________________________________________________________________
/**
 * Constructor.
 */
OutOfRangeError::OutOfRangeError()
{
	;
}

/**
 * Constructor.
 */
OutOfRangeError::OutOfRangeError(const std::string & details)
{
	std::ostringstream oss;
	oss << "mfm::OutOfRangeError: " << details;
	message = oss.str();
}

/**
 * Destructor.
 */
OutOfRangeError::~OutOfRangeError() throw()
{
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string OutOfRangeError::formatMessage() const
{
	std::ostringstream message;
	message << "mfm::OutOfRangeError";
	return message.str();
}
//______________________________________________________________________
/**
 * Constructor.
 */
FieldNotFound::FieldNotFound(std::string const & _fieldName)
	: fieldName_(_fieldName)
{
	;
}

/**
 * Destructor.
 */
FieldNotFound::~FieldNotFound() throw()
{
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string FieldNotFound::formatMessage() const
{
	std::ostringstream message;
	message << "Could not find field " << fieldName_ << std::endl;
	return message.str();
}
//______________________________________________________________________
/**
 * Constructor.
 */
BitFieldNotFound::BitFieldNotFound(std::string const & _fieldName, std::string const & _bitFieldName)
	: fieldName_(_fieldName), bitFieldName_(_bitFieldName)
{
	;
}

/**
 * Destructor.
 */
BitFieldNotFound::~BitFieldNotFound() throw()
{
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string BitFieldNotFound::formatMessage() const
{
	std::ostringstream message;
	message << "Could not find bit field " << bitFieldName_ << " within field '" << fieldName_ << "'" << std::endl;
	return message.str();
}
//______________________________________________________________________
/**
 * Constructor.
 */
FormatNotFound::FormatNotFound(uint16_t type)
	: type_(type)
{
	;
}

/**
 * Destructor.
 */
FormatNotFound::~FormatNotFound() throw()
{
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string FormatNotFound::formatMessage() const
{
	std::ostringstream message;
	message << "No frame format found with type '" << type_ << "'" << std::endl;
	return message.str();
}
//______________________________________________________________________
/**
 * Constructor.
 */
FormatRevisionNotFound::FormatRevisionNotFound(uint16_t type, uint8_t revision)
	: FormatNotFound(type), revision_(revision)
{
	;
}

/**
 * Destructor.
 */
FormatRevisionNotFound::~FormatRevisionNotFound() throw()
{
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string FormatRevisionNotFound::formatMessage() const
{
	std::ostringstream message;
	message << "No revision '" << (short) revision_ << "' found for format with type '" << type_ << "'" << std::endl;
	return message.str();
}
//______________________________________________________________________
} /* namespace mfm */
