// =============================================================================
// IncludeHandler.cpp
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

#include "IncludeHandler.h"
#include "../Object.h"

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::IncludeHandler class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
IncludeHandler::IncludeHandler()
{
}

/**
 * Destructor.
 */
IncludeHandler::~IncludeHandler()
{
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::NoOpIncludeHandler class methods body
///////////////////////////////////////////////////////////////////////////////

std::auto_ptr<NoOpIncludeHandler> NoOpIncludeHandler::instance;

/**
 * Constructor.
 */
NoOpIncludeHandler::NoOpIncludeHandler()
{
}

/**
 * Destructor.
 */
NoOpIncludeHandler::~NoOpIncludeHandler()
{
}

/**
 * No-operation callback method to ignore includes.
 */
void NoOpIncludeHandler::parseInclude(Object&)
{
}

/**
 * No-operation callback method to ingore includes.
 * @return Returns always false, which indicates to go on serializing the tree
 * as if there was no include.
 */
bool NoOpIncludeHandler::formatInclude(const Object&)
{
	return false;
}

/**
 * Singleton getter.
 * @return Returns the instance of the singleton.
 */
NoOpIncludeHandler& NoOpIncludeHandler::getInstance()
{
	if(instance.get() == 0)
	{
		instance.reset(new NoOpIncludeHandler());
	}
	return *instance;
}

}

}
