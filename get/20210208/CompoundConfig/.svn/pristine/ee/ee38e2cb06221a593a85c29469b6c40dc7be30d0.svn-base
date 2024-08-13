// =============================================================================
// FileStorageFactory.cpp
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

#include "FileStorageFactory.h"
#include "FileStorage.h"
#include "Url.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::FileStorageFactory class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * (see inherited doc)
 */
std::auto_ptr<Storage> FileStorageFactory::createStorage(const Url& url) const
{
	std::auto_ptr<Storage> storage;
	if(url.isFile())
	{
		storage.reset(new FileStorage);
		// we don't define the path in the constructor because we don't want to
		// open the file yet (besides, we don't know the open mode either).
		static_cast<FileStorage&>(*storage).setFileName(url.toPath());
	}
	return storage;
}

}
}
