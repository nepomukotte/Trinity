// =============================================================================
// TextFormat.cpp
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

#include "TextFormat.h"
#include "FileStorage.h"
#include "Url.h"

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::TextFormat class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
TextFormat::TextFormat()
{
}

/**
 * Destructor.
 */
TextFormat::~TextFormat()
{
}

std::string TextFormat::getConfigName(Storage& storage)
{
	FileStorage* file = dynamic_cast<FileStorage*>(&storage);
	if(file != 0)
	{
		const std::string& path = file->getFileName();
		size_t pos = path.find_last_of(Url::PATH_SEPARATOR);
		if(pos != std::string::npos)
		{
			return path.substr(pos + 1);	// return the substring of the file name.
		}
		else // probably only a file name in the current directory.
		{
			return path;
		}
	}
	return "unnamed";

}

}
}
