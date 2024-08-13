// =============================================================================
// FileStorage.cpp
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

#include "FileStorage.h"
#include "Url.h"
#include <sstream>
#include <algorithm>
#include <cassert>

using std::ios;

#undef min

namespace CCfg
{

namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::FileStorage class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Convert a FileStorage open mode to a STL open mode.
 * @param mode the FileStorage open mode to convert.
 * @return Returns the converted STL open mode.
 */
ios::openmode FileStorage::convertOpenMode(OpenMode mode)
{
	std::ios::openmode stlMode = ios::in;
	if(mode & Read)		stlMode |= ios::in;
	else				stlMode ^= ios::in;
	if(mode & Write)	stlMode |= ios::out;
	if(mode & Append)	stlMode |= ios::app;
	if(mode & Truncate)	stlMode |= ios::trunc;
	if(mode & Binary)	stlMode |= ios::binary;

	return stlMode;
}

/**
 * Default constructor.
 */
FileStorage::FileStorage() : fileName(), url(), openMode(), file()
{
}

/**
 * Full constructor. Opens a file.
 * @param fileName the name of the file to open.
 * @param mode the open mode.
 * @throws FileStorage::FileOpenError if the file cannot be opened.
 */
FileStorage::FileStorage(const std::string& fileName, OpenMode mode)
	: fileName(fileName), url(Url::PathNormalizer::escape(fileName)),
	openMode(mode), file(fileName.c_str(), convertOpenMode(mode))
{
	if(file.fail())
	{
		throw FileOpenError(fileName, openMode);
	}
}

/**
 * Destructor.
 */
FileStorage::~FileStorage()
{
	// fstream::close is automatically called. No need to call FileStorage::close.
}

/**
 * Getter of the "fileName" property.
 * @return Returns the name of the currently opened file.
 */
const std::string& FileStorage::getFileName() const
{
	return fileName;
}

/**
 * Getter of the "openMode" property.
 * @return Returns the open mode of the currently opened file.
 */
FileStorage::OpenMode FileStorage::getOpenMode() const
{
	return openMode;
}

/**
 * Setter of the "fileName" property.
 * @param fileName the name of the file to open in the next call to open().
 * If a file is already opened, the property won't change.
 */
void FileStorage::setFileName(const std::string& fileName)
{
	if(not isOpened())
	{
		this->fileName = fileName;
		this->url = Url::PathNormalizer::escape(fileName);
	}
}

/**
 * Setter of the "openMode" property.
 * @param mode the open mode that the next call to open() will use.
 * If a file is already opened, the property won't change.
 */
void FileStorage::setOpenMode(OpenMode mode)
{
	if(not isOpened())
	{
		openMode = mode;
	}
}

/**
 * Opens a file using fileName and openMode internal properties.
 * @throws FileStorage::FileOpenError if the file cannot be opened.
 */
void FileStorage::open()
{
	open(fileName, openMode);
}

/**
 * Opens the file using fileName property in read or write mode (or both).
 * The internal open mode will only keep its flags "Binary".
 * @param write enable write mode. if the file exists, it will be overwritten.
 */
void FileStorage::open(bool write)
{
	openMode &= Binary; // clears open mode, except Binary flag
	if(write)
	{
		openMode |= Write | Truncate;
	}
	else
	{
		openMode |= Read;
	}

	open(fileName, openMode);
}

/**
 * Opens a file.
 * @param fileName the name of the file to open.
 * @param mode the open mode.
 * @throws FileStorage::FileOpenError if the file cannot be opened.
 */
void FileStorage::open(const std::string& fileName, OpenMode mode)
{
	if(file.is_open())
	{
		file.close();
	}
	setFileName(fileName);
	setOpenMode(mode);
	file.open(fileName.c_str(), convertOpenMode(openMode));
	if(file.fail())
	{
		throw FileOpenError(fileName, openMode);
	}
}

/**
 * Closes the file.
 * If no file is opened, it has no effect.
 */
void FileStorage::close()
{
	file.close();
}

/**
 * Indicates whether the file is opened or not.
 */
bool FileStorage::isOpened() const
{
	return const_cast<std::fstream*>(&file)->is_open();
}

/*
 *
 */
const std::string& FileStorage::getUrl() const
{
	return url;
}

/**
 * Reads the file completely.
 * @param data the buffer for read data.
 */
void FileStorage::readAll(std::string& data)
{
	char buffer[128];
	do
	{
		file.read(buffer, sizeof buffer);
		data.append(buffer, file.gcount());
	}
	while(not file.eof());
}

/**
 * Reads some data in the file.
 * @param count the number of bytes to read.
 * @param data the buffer for read data.
 */
void FileStorage::read(size_t count, std::string& data)
{
	size_t readCount = 0;
	char buffer[128];
	do
	{
		file.read(buffer, std::min(sizeof buffer, count - readCount));
		data.append(buffer, file.gcount());
		readCount += file.gcount();
	}
	while(readCount < count);
}

/**
 * Writes a string into the file.
 * @param data the string to write.
 */
void FileStorage::write(const std::string& data)
{
	file.write(data.data(), data.size());
}

/**
 * Writes some binary data into the buffer.
 */
void FileStorage::write(const void* data, unsigned sizeBytes)
{
	file.write(static_cast<const char*>(data), sizeBytes);
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::FileStorage::FileOpenError class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param fileName the name of the file which cannot be opened.
 * @param openMode the mode in which we tried to open the file.
 */
FileStorage::FileOpenError::FileOpenError(const std::string& fileName, OpenMode openMode)
	: fileName(fileName), openMode(openMode)
{
}

/**
 * Copy constructor: initializes a new exception instance with a copy of the
 * informations in the exception in parameter.
 * @param value the exception to copy.
 */
FileStorage::FileOpenError::FileOpenError(const FileOpenError& value)
	: OpenError(value), fileName(value.fileName), openMode(value.openMode)
{

}

/**
 * Destructor.
 */
FileStorage::FileOpenError::~FileOpenError() throw()
{
}

/**
 * Assignment operator.
 * @param value the exception to copy
 */
FileStorage::FileOpenError& FileStorage::FileOpenError::operator=(const FileOpenError& value)
{
	OpenError::operator =(value);
	fileName = value.fileName;
	openMode = value.openMode;
	return *this;
}

/**
 * Getter of the "fileName" property.
 * @return Returns the name of the file which cannot be opened.
 */
const std::string& FileStorage::FileOpenError::getFileName() const
{
	return fileName;
}

/**
 * Getter of the "openMode" property.
 * @return Returns the mode in which we tried to open the file.
 */
FileStorage::OpenMode FileStorage::FileOpenError::getOpenMode() const
{
	return openMode;
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string FileStorage::FileOpenError::formatMessage() const
{
	std::ostringstream message;
	message << "Cannot open the file \"" << fileName
		<< "\" using mode: " << openMode;
	return message.str();
}


}

}
