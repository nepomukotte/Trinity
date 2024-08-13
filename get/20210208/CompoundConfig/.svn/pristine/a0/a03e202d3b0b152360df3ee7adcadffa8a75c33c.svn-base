// =============================================================================
// StorageFactory.h
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

#ifndef CCfg_Io_StorageFactory_h_INCLUDED
#define CCfg_Io_StorageFactory_h_INCLUDED

#include "AbstractStorageFactory.h"
#include <boost/ptr_container/ptr_vector.hpp>

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class is a singleton that creates storages from URLs using a list of
 * registered factories.
 * The file storage is automatically registered.
 *
 * @author fchateau
 */
class StorageFactory : public AbstractStorageFactory
{
	boost::ptr_vector<AbstractStorageFactory> factories;

public:
	~StorageFactory() {}
	std::auto_ptr<Storage> createStorage(const Url& url) const;
	void registerStorage(std::auto_ptr<AbstractStorageFactory> factory);

	static StorageFactory& instance();

private:
	StorageFactory();
	StorageFactory(const StorageFactory& rvalue);
	StorageFactory& operator=(const StorageFactory& rvalue);
};

}
}

#endif // CCfg_Io_StorageFactory_h_INCLUDED
