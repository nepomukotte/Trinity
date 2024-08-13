// =============================================================================
// StreamStorage.h
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

#ifndef CCfg_Io_StreamStorage_h_INCLUDED
#define CCfg_Io_StreamStorage_h_INCLUDED

#include "Storage.h"
#include <string>

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * An intermediary interface to distinguish storage families.
 *
 * @author fchateau
 */
class StreamStorage : public Storage
{
public:
	StreamStorage();
	virtual ~StreamStorage();

	virtual void open() = 0;
	virtual void open(bool write) = 0;
	virtual void close() = 0;
	virtual bool isOpened() const = 0;

	virtual void readAll(std::string& data) = 0;
	virtual void read(size_t count, std::string& data) = 0;
	virtual void write(const std::string& data) = 0;
	virtual void write(const void* data, unsigned sizeBytes) = 0;
};

}

}

#endif // CCfg_Io_StreamStorage_h_INCLUDED
