// =============================================================================
// UnhandledUrl.cpp
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

#include "UnhandledUrl.h"
#include <sstream>

namespace CCfg
{
namespace Io
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::UnhandledUrl class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param url the unhandled URL.
 * @param product the name of the product whose creation failed.
 */
UnhandledUrl::UnhandledUrl(const std::string& url, const std::string& product)
	: url(url), product(product)
{
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string UnhandledUrl::formatMessage() const
{
	std::ostringstream message;
	message << "Cannot build a " << product << " from the URL: \"" << url << "\".";
	return message.str();
}

}
}
