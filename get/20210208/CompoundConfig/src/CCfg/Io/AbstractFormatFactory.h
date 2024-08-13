// =============================================================================
// AbstractFormatFactory.h
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

#ifndef CCfg_Io_AbstractFormatFactory_h_INCLUDED
#define CCfg_Io_AbstractFormatFactory_h_INCLUDED

#include <memory>

namespace CCfg
{
namespace Io
{
// forward declarations
class Format;
class Url;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is an abstract factory that builds configuration formats from Urls.
 * To integrate a new format in the framework, a subclass of this factory must
 * be created.
 *
 * @author fchateau
 */
class AbstractFormatFactory
{
public:
	virtual ~AbstractFormatFactory() {}
	virtual std::auto_ptr<Format> createFormat(const Url& url) const = 0;
};

}
}

#endif // CCfg_Io_AbstractFormatFactory_h_INCLUDED
