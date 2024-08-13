// =============================================================================
// StorageFactory.cpp
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

#include "StorageFactory.h"
#include "FileStorageFactory.h"
#include "UnhandledUrl.h"
#include "Storage.h"
#include "Url.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::StorageFactory class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Private constructor.
 * Automatically registers the file storage factory.
 */
StorageFactory::StorageFactory()
{
	registerStorage(std::auto_ptr<AbstractStorageFactory>(new FileStorageFactory));
}

/*
 * (see inherited doc)
 */
std::auto_ptr<Storage> StorageFactory::createStorage(const Url& url) const
{
	std::auto_ptr<Storage> storage;
	for(boost::ptr_vector<AbstractStorageFactory>::const_iterator curFactory = factories.begin(),
		endFactory = factories.end();
		curFactory != endFactory and storage.get() == 0;
		++curFactory)
	{
		std::auto_ptr<Storage> storage2 = curFactory->createStorage(url);
		storage = storage2; // hack for vxWorks
	}
	if(storage.get() == 0)
	{
		throw UnhandledUrl(url.format(), "Storage");
	}
	return storage;
}

/**
 * Registers a storage factory to enable the recognition of a storage in the whole
 * application.
 */
void StorageFactory::registerStorage(std::auto_ptr<AbstractStorageFactory> factory)
{
	factories.push_back(factory.release());
}

/**
 * Returns the singleton instance.
 */
StorageFactory& StorageFactory::instance()
{
	static StorageFactory inst;
	return inst;
}

}
}
