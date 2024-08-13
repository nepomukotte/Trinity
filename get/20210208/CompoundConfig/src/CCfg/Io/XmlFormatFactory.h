// =============================================================================
// XmlFormatFactory.h
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

#ifndef CCfg_Io_XmlFormatFactory_h_INCLUDED
#define CCfg_Io_XmlFormatFactory_h_INCLUDED

#include "AbstractFormatFactory.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class creates XmlFormat objects from URLs.
 *
 * @author fchateau
 */
class XmlFormatFactory : public AbstractFormatFactory
{
	static const int XML_EXT_COUNT = 2;
	static const char* const XML_EXTENSIONS[XML_EXT_COUNT];

public:
	~XmlFormatFactory() {}
	std::auto_ptr<Format> createFormat(const Url& url) const;
};

}
}

#endif // CCfg_Io_XmlFormatFactory_h_INCLUDED
