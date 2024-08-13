// =============================================================================
// FileStorageFactory.h
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

#ifndef CCfg_Io_FileStorageFactory_h_INCLUDED
#define CCfg_Io_FileStorageFactory_h_INCLUDED

#include "AbstractStorageFactory.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class is a factory that produces Storage objects to store configurations in files.
 *
 * @author fchateau
 */
class FileStorageFactory : public AbstractStorageFactory
{
public:
	~FileStorageFactory() {}
	std::auto_ptr<Storage> createStorage(const Url& url) const;
};

}
}

#endif // CCfg_Io_FileStorageFactory_h_INCLUDED
