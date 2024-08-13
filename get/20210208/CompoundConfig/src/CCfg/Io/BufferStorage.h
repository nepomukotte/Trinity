// =============================================================================
// BufferStorage.h
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

#ifndef CCfg_Io_BufferStorage_h_INCLUDED
#define CCfg_Io_BufferStorage_h_INCLUDED

#include "StreamStorage.h"
#include "CCfg/Object.h"
#include "CCfg/Io/IncludeHandler.h"
#include <sstream>

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * Abstraction layer for memory buffer storage.
 *
 * @author fchateau
 */
class BufferStorage : public StreamStorage
{
	std::stringstream buffer;

public:
	// construction/destruction
	BufferStorage();
	virtual ~BufferStorage();

	// opening
	void open();
	void open(bool write);
	void close();
	bool isOpened() const;

	// Url
	const std::string& getUrl() const;

	// input/outputs
	void readAll(std::string& buffer);
	void read(size_t count, std::string& buffer);
	void write(const std::string& data);
	void write(const void* data, unsigned sizeBytes);

	template<class F>
	static std::string encode(Object* ccfgObj);

	template<class F>
	static std::auto_ptr<Object> decode(const std::string& ccfgStr);
};

/**
 * Serializer fonction for CConfig trees.
 * @tparam F Format class
 * @param ccfgObj The CompoundCOnfig object to be serialized.
 * @return A std::string containing the serialized data.
 */
template<class F>
std::string BufferStorage::encode(Object* ccfgObj)
{
	if (ccfgObj == 0) return "";
	CCfg::Io::BufferStorage bufStore;
	F fmt;
	fmt.format(ccfgObj, bufStore, CCfg::Io::NoOpIncludeHandler::getInstance());
	return bufStore.buffer.str();
}

/**
 * De-serializer fonction for CConfig trees.
 * @tparam F Format class
 * @param ccfgStr The std::string containing the serialized CompoundConfig to be de-serialized.
 * @return The de-serialized CCfg::Object auto_ptr.
 */
template<class F>
std::auto_ptr<Object> BufferStorage::decode(const std::string& ccfgStr)
{
	CCfg::Io::BufferStorage bufStore;
	bufStore.write(ccfgStr);
	F fmt;
	return fmt.parse(bufStore, CCfg::Io::NoOpIncludeHandler::getInstance());
}

/**
 * Helper serializer template function. @see BufferStorage::encode
 */
template<class F>
inline std::string encode(Object* ccfgObj)
{
	return BufferStorage::encode<F>(ccfgObj);
}

/**
 * Helper de-serializer template function. @see BufferStorage::decode
 */
template<class F>
std::auto_ptr<Object> decode(const std::string& ccfgStr)
{
	return BufferStorage::decode<F>(ccfgStr);
}

}
}

#endif // CCfg_Io_BufferStorage_h_INCLUDED
