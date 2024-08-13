// =============================================================================
// BufferStorage.cpp
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

#include "BufferStorage.h"
#include <cassert>

#undef min

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::BufferStorage class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
BufferStorage::BufferStorage()
	: buffer()
{
}

/**
 * Destructor.
 */
BufferStorage::~BufferStorage()
{
}

/**
 * Resets the buffer.
 */
void BufferStorage::open()
{
	//buffer.str("");
}

/**
 * Resets the buffer.
 */
void BufferStorage::open(bool)
{
	open();
}

/**
 * Resets the buffer.
 */
void BufferStorage::close()
{
	//buffer.str("");
}

/**
 * Always returns true.
 */
bool BufferStorage::isOpened() const
{
	return true;
}

/**
 * Returns an empty string because a buffer is a transient resource.
 */
const std::string& BufferStorage::getUrl() const
{
	static std::string dummyUrl;
	return dummyUrl;
}


/**
 * Reads the buffer completely.
 * @param data the placeholder for read data.
 */
void BufferStorage::readAll(std::string& data)
{
	data = buffer.str();
}

/**
 * Reads some data in the buffer.
 * @param count the number of bytes to read.
 * @param data the placeholder for read data.
 */
void BufferStorage::read(size_t count, std::string& data)
{
	size_t readCount = 0;
	char buf[128];
	do
	{
		buffer.read(buf, std::min(sizeof buf, count - readCount));
		data.append(buf, buffer.gcount());
		readCount += buffer.gcount();
	}
	while(readCount < count);
}

/**
 * Writes a string into the buffer.
 * @param data the data to write.
 */
void BufferStorage::write(const std::string& data)
{
	buffer.write(data.data(), data.size());
}

/**
 * Writes some binary data into the buffer.
 */
void BufferStorage::write(const void* data, unsigned sizeBytes)
{
	buffer.write(static_cast<const char*>(data), sizeBytes);
}

}
}
