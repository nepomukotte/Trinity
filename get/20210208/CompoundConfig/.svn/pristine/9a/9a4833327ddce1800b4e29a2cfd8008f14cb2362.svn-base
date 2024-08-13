// =============================================================================
// Storage.h
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

#ifndef CCfg_Io_Storage_h_INCLUDED
#define CCfg_Io_Storage_h_INCLUDED

#include "../Exception.h"

namespace CCfg
{
namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * An interface that each storage facility must realize.
 *
 * @author fchateau
 */
class Storage
{
public:
	// inner classes
	class OpenError;

	// construction/destruction
	Storage();
	virtual ~Storage();

	// opening
	virtual void open() = 0;
	virtual void close() = 0;
	virtual bool isOpened() const = 0;

	// Url
	virtual const std::string& getUrl() const = 0;

private:
	// non-copyable
	Storage(const Storage&);
	Storage& operator=(const Storage&);
};


/**
 * Signals an error while opening a storage.
 *
 * @author fchateau
 */
class Storage::OpenError : public Exception
{
public:
	// construction/destruction
	OpenError();
	OpenError(const std::string& message);
	OpenError(const OpenError& value);
	virtual ~OpenError() throw();
	OpenError& operator=(const OpenError& value);
};

}

}

#endif // CCfg_Io_Storage_h_INCLUDED
