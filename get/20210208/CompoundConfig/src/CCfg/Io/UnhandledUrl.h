// =============================================================================
// UnhandledUrl.h
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

#ifndef CCfg_Io_UnhandledUrl_h_INCLUDED
#define CCfg_Io_UnhandledUrl_h_INCLUDED

#include "../Exception.h"
#include <string>

namespace CCfg
{
namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This exception reports that a factory was unable to build a product for the
 * specified URL.
 */
class UnhandledUrl : public Exception
{
	std::string url;		///< the unhandled URL.
	std::string product;	///< the name of the product whose creation failed.

public:
	// construction/destruction
	UnhandledUrl(const std::string& url, const std::string& product);
	virtual ~UnhandledUrl() throw() {}

	// getters
	const std::string& getUrl() const { return url; }
	const std::string& getProduct() const { return product; }

private:
	std::string formatMessage() const;
};

}
}

#endif // CCfg_Io_UnhandledUrl_h_INCLUDED
