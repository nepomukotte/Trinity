// =============================================================================
// UnsupportedOperation.cpp
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

#include "UnsupportedOperation.h"

#include <sstream>

using std::string;
using std::ostringstream;
using std::type_info;

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// UnsupportedOperation class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor: initializes informations about the exception.
 * @param name \copydoc #name
 * @param clazz the type_info describing the class that reports the error.
 * @param operationName the name of the unsupported operation.
 * @param comment an additionnal comment about the error.
 */
UnsupportedOperation::UnsupportedOperation(
	const std::string& name,
	const std::type_info& clazz,
	const std::string& operationName,
	const std::string& comment)
	: name(name), clazz(clazz), operationName(operationName), comment(comment)
{
}

/**
 * Copy constructor: initializes a new exception instance with a copy of the
 * informations in the exception in parameter.
 * @param value the exception to copy.
 */
UnsupportedOperation::UnsupportedOperation(const UnsupportedOperation& value) :
		Exception(value),
		name(value.name),
		clazz(value.clazz),
		operationName(value.operationName),
		comment(value.comment)
{
}

/**
 * Destructor.
 */
UnsupportedOperation::~UnsupportedOperation() throw()
{
}

std::string UnsupportedOperation::formatMessage() const
{
	ostringstream message;
	message << "The operation: \"" << operationName << "\" is unsupported "
			<< "by the object \'" << name << '\''
			<< " of type " << clazz.name() << ". "
			<< comment;
	return message.str();
}

}

