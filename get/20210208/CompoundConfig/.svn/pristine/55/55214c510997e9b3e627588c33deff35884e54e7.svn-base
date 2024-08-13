// =============================================================================
// Storage.cpp
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

#include "Storage.h"
#include <assert.h>

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Storage class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
Storage::Storage()
{
}

/**
 * Destructor.
 */
Storage::~Storage()
{
}

/**
 * @fn Storage::open()
 * Opens the storage.
 * @throws Storage::OpenError if the file cannot be opened.
 */

/**
 * @fn Storage::close()
 * Closes the storage.
 */

/**
 * @fn Storage::isOpened() const
 * Returns true if the storage is currently opened.
 */

/**
 * @fn Storage::getUrl() const
 * Returns the URL to the storage.
 */

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Storage::OpenError class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
Storage::OpenError::OpenError()
	: Exception()
{
}

/**
 * Constructor from an error message.
 */
Storage::OpenError::OpenError(const std::string& message)
	: Exception(message)
{
}

/**
 * Copy constructor.
 */
Storage::OpenError::OpenError(const OpenError& value)
	: Exception(value)
{
}

/**
 * Destructor.
 */
Storage::OpenError::~OpenError() throw()
{
}

/**
 * Assignment operator.
 */
Storage::OpenError& Storage::OpenError::operator=(const OpenError& value)
{
	Exception::operator=(value);
	return *this;
}

}

}
