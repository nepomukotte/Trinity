// =============================================================================
// XmlFormatFactory.cpp
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

#include "XmlFormatFactory.h"
#include "XmlFormat.h"
#include "Url.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::XmlFormatFactory class methods body
////////////////////////////////////////////////////////////////////////////////

const int XmlFormatFactory::XML_EXT_COUNT;
const char* const XmlFormatFactory::XML_EXTENSIONS[XML_EXT_COUNT] = {"xml", "xcfg"};

/*
 * (see inherited doc)
 */
std::auto_ptr<Format> XmlFormatFactory::createFormat(const Url& url) const
{
	std::auto_ptr<Format> format;
	if(url.isFile())
	{
		std::string extension = url.getExtension();

		for(int i = 0; i < XML_EXT_COUNT; i++)
		{
			if(extension ==	XML_EXTENSIONS[i])
			{
				format.reset(new XmlFormat);
			}
		}
	}
	return format;
}

}
}
