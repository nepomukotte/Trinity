// =============================================================================
// Format.cpp
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

#include "Format.h"
#include "Storage.h"
#include "IncludeHandler.h"
#include "../Object.h"


namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Format class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
Format::Format()
{
}

/**
 * Destructor.
 */
Format::~Format()
{
}

/**
 * @fn std::auto_ptr<CCfg::Object> Format::parse(Storage& storage, IncludeHandler& handler)
 * Parses the storage data to produce a configuration tree.
 * @param storage the storage where is the data.
 * @param handler a callback to react to include directives.
 * @return Returns the created configuration tree.
 * @throws Storage::OpenError if there was an error opening the storage.
 */

/**
 * @fn void Format::format(const Object *config, Storage& storage, IncludeHandler& handler)
 * Writes a configuration tree in a storage.
 * @param config the configuration tree.
 * @param storage the storage where data is written.
 * @param handler a callback to react to include directives.
 * @throws Storage::OpenError if there was an error opening the storage.
 */

/**
 * @fn std::string Format::getConfigName(Storage& storage)
 * Returns the name of the configuration in a storage.
 * @param storage the storage of the configuration.
 * @return Returns the name of the configuration.
 */

/**
 * Convenience method to parse a configuration while ignoring includes.
 */
std::auto_ptr<CCfg::Object> Format::parse(Storage& storage)
{
	return parse(storage, NoOpIncludeHandler::getInstance());
}

/**
 * Convenience method to format a configuration while ignoring includes.
 */
void Format::format(const Object *config, Storage& storage)
{
	format(config, storage, NoOpIncludeHandler::getInstance());
}


///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Format::InvalidFormat class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
Format::InvalidFormat::InvalidFormat()
{
}

/**
 * Destructor.
 */
Format::InvalidFormat::~InvalidFormat() throw()
{
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Format::SyntaxError class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
Format::SyntaxError::SyntaxError()
{
}

/**
 * Destructor.
 */
Format::SyntaxError::~SyntaxError() throw()
{
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Format::SemanticError class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
Format::SemanticError::SemanticError()
{
}

/**
 * Destructor.
 */
Format::SemanticError::~SemanticError() throw()
{
}

}

}
