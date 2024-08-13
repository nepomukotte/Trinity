// =============================================================================
// AbstractStorageFactory.h
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

#ifndef CCfg_Io_AbstractStorageFactory_h_INCLUDED
#define CCfg_Io_AbstractStorageFactory_h_INCLUDED

#include <memory>

namespace CCfg
{
namespace Io
{
// forward declarations
class Storage;
class Url;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is an abstract factory that builds storage systems for configurations from URLs.
 * To integrate a new storage in the framework, a subclass of this factory must
 * be created.
 *
 * @author fchateau
 */
class AbstractStorageFactory
{
public:
	virtual ~AbstractStorageFactory() {}
	virtual std::auto_ptr<Storage> createStorage(const Url& url) const = 0;
};

}
}

#endif // CCfg_Io_AbstractStorageFactory_h_INCLUDED
