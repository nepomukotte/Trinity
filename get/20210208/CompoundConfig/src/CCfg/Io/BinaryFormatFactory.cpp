// =============================================================================
// BinaryFormatFactory.cpp
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Patrick Sizun (patrick.sizun@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "BinaryFormatFactory.h"
#include "BinaryFormat.h"
#include "Url.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::BinaryFormatFactory class methods body
////////////////////////////////////////////////////////////////////////////////

const int BinaryFormatFactory::BINARY_EXT_COUNT;
const char* const BinaryFormatFactory::BINARY_EXTENSIONS[BINARY_EXT_COUNT] = {"bin", "dat", "bcfg"};

/*
 * (see inherited doc)
 */
std::auto_ptr<Format> BinaryFormatFactory::createFormat(const Url& url) const
{
	std::auto_ptr<Format> format;
	if(url.isFile())
	{
		std::string extension = url.getExtension();

		for(int i = 0; i < BINARY_EXT_COUNT; i++)
		{
			if(extension ==	BINARY_EXTENSIONS[i])
			{
				format.reset(new BinaryFormat);
			}
		}
	}
	return format;
}

}
}
