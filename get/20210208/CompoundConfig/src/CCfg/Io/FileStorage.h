// =============================================================================
// FileStorage.h
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

#ifndef CCfg_Io_FileStorage_h_INCLUDED
#define CCfg_Io_FileStorage_h_INCLUDED

#include "StreamStorage.h"
#include <string>
#include <fstream>

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * Abstraction layer for file storage.
 *
 * @author fchateau
 */
class FileStorage : public StreamStorage
{
public:
	// inner classes
	class FileOpenError;

	// enums/typedefs
	enum OpenFlags {Read = 1, Write = 2, Append = 4, Truncate = 8, Binary = 16};
	typedef unsigned OpenMode;

	// construction/destruction
	FileStorage();
	FileStorage(const std::string& fileName, OpenMode mode);
	virtual ~FileStorage();

	// getters/setters
	void setFileName(const std::string& fileName);
	void setOpenMode(OpenMode mode);
	const std::string& getFileName() const;
	OpenMode getOpenMode() const;

	// opening
	void open();
	void open(bool write);
	void open(const std::string& fileName, OpenMode mode);
	void close();
	bool isOpened() const;

	// Url
	const std::string& getUrl() const;

	// input/outputs
	void readAll(std::string& data);
	void read(size_t count, std::string& data);
	void write(const std::string& data);
	void write(const void* data, unsigned sizeBytes);

private:
	static std::ios::openmode convertOpenMode(OpenMode mode);

	// attributes
	std::string fileName;	///< the name of the opened file.
	std::string url;		///< the url of the opened file.
	OpenMode openMode;		///< the open mode of the opened file.
	std::fstream file;		///< the file stream.
};

/**
 * Signals an error while opening a file.
 *
 * @author fchateau
 */
class FileStorage::FileOpenError : public Storage::OpenError
{
	std::string fileName;	///< the name of the file which cannot be opened.
	OpenMode openMode;		///< the mode in which we tried to open the file.
public:
	// construction/destruction
	FileOpenError(const std::string& fileName, OpenMode openMode);
	FileOpenError(const FileOpenError& value);
	virtual ~FileOpenError() throw();
	FileOpenError& operator=(const FileOpenError& value);

	// getters
	const std::string& getFileName() const;
	OpenMode getOpenMode() const;

private:
	std::string formatMessage() const;
};

}

}

#endif // CCfg_Io_FileStorage_h_INCLUDED
