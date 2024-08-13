// =============================================================================
// Document.cpp
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

#include "Document.h"
#include "FormatFactory.h"
#include "StorageFactory.h"
#include "UnhandledUrl.h"
#include "Storage.h"
#include "Format.h"
#include "Url.h"
#include "../Object.h"
#include <sstream>
#include <cassert>

#ifdef VXWORKS
#include <ioLib.h>
#else
#include <unistd.h>
#endif


namespace CCfg
{
namespace Io
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Document class methods body
///////////////////////////////////////////////////////////////////////////////

#define RETHROW_AS_IO_ERROR(type, url)\
	catch(Storage::OpenError& e)\
	{\
		throw IoError(type, url, e.getMessage());\
	}\
	catch(Format::InvalidFormat& e)\
	{\
		throw IoError(type, url, e.getMessage());\
	}\
	catch(Url::ParseError& e)\
	{\
		throw IoError(type, url, e.getMessage());\
	}\
	catch(UnhandledUrl& e)\
	{\
		throw IoError(type, url, e.getMessage());\
	}\
	catch(Document::RecursiveInclude& e)\
	{\
		throw IoError(type, url, e.getMessage());\
	}

/**
 * Default constructor. Creates an uninitialized document.
 */
Document::Document()
	: config(new Object), name("new config")
{
}

/**
 * Implicit constructor. Creates a document and loads the specified resource.
 * @param srcUrl the location of the resource to load.
 */
Document::Document(const Url& srcUrl)
	: config(), name()
{
	loadFromUrl(srcUrl);
}

Document::Document(const std::string& path)
	: config(), name()
{
	loadFromFile(path);
}

/**
 * Destructor. Destroy the owned configuration.
 */
Document::~Document()
{
}

/**
 * Loads a configuration from a file.
 * @throws Document::IoError reports that the document cannot be loaded.
 */
Object& Document::loadFromFile(const std::string& path, bool ignoreIncludes)
{
	return loadFromUrl(Url::fromPath(path), ignoreIncludes);
}

/**
 * Saves a configuration to a file
 * @throws Document::IoError reports that the document cannot be saved.
 */
void Document::saveToFile(const std::string& path, bool flatten)
{
	saveToUrl(Url::fromPath(path), flatten);
}

/**
 * Loads a configuration from a storage in the document.
 * @param srcUrl the URL of the storage.
 * @param ignoreIncludes indicates whether include directive should be ignored.
 * Default is false: all includes are done.
 * @return Returns the configuration root object.
 * @throws Document::IoError reports that the document cannot be loaded.
 */
Object& Document::loadFromUrl(const Url& srcUrl, bool ignoreIncludes)
{
	try
	{
		std::auto_ptr<Format> format = FormatFactory::instance().createFormat(srcUrl);
		std::auto_ptr<Storage> storage = StorageFactory::instance().createStorage(srcUrl);

		return load(*storage, *format, ignoreIncludes);
	}
	RETHROW_AS_IO_ERROR(IoError::Load, srcUrl.format());
}

/**
 * Saves the configuration document in a unique storage (ignoring includes).
 * @param destUrl the URL of the storage.
 * @param flatten If true, all included files are merged into a flat configuration.
 * @throws Document::IoError reports that the document cannot be saved.
 */
void Document::saveToUrl(const Url& destUrl, bool flatten)
{
	try
	{
		std::auto_ptr<Format> format = FormatFactory::instance().createFormat(destUrl);
		std::auto_ptr<Storage> storage = StorageFactory::instance().createStorage(destUrl);

		save(*storage, *format, flatten);
	}
	RETHROW_AS_IO_ERROR(IoError::Save, destUrl.format());
}

/**
 * Loads a configuration from a storage/format in the document.
 * @param srcStorage the storage where the configuration should be read.
 * @param format the configuration decoding format.
 * @param ignoreIncludes indicates whether include directive should be ignored.
 * Default is false: all includes are done.
 * @throws Document::IoError reports that the document cannot be loaded.
 */
Object& Document::load(Storage& srcStorage, Format& format, bool ignoreIncludes)
{
	try
	{
		IncludeHandler& handler = ignoreIncludes ?
			NoOpIncludeHandler::getInstance() : static_cast<IncludeHandler&>(*this);

		Url url(srcStorage.getUrl());
		includeStack.clear();
		includeStack.push_back(&url);

		std::auto_ptr<Object> temp = load(srcStorage, format, handler);
		config = temp;	// hack for vxworks which doesn't support assignment from return!!

		includeStack.pop_back();

		if(config.get() == 0)
		{
			throw IoError(IoError::Load, srcStorage.getUrl(), "empty configuration!");
		}

		setName(format.getConfigName(srcStorage));
		setUrl(srcStorage.getUrl());
		return *config;
	}
	RETHROW_AS_IO_ERROR(IoError::Load, srcStorage.getUrl());
}

/**
 * Saves the configuration document in a unique storage (ignoring includes).
 * @param destStorage the URL of the storage.
 * @param format the configuration encoding format.
 * @param flatten If true, all included files are merged into a flat configuration.
 * @throws Document::IoError reports that the document cannot be saved.
 */
void Document::save(Storage& destStorage, Format& format, bool flatten)
{
	try
	{
		IncludeHandler& handler = flatten ?
			NoOpIncludeHandler::getInstance() : static_cast<IncludeHandler&>(*this);

		Url url(destStorage.getUrl());
		includeStack.clear();
		includeStack.push_back(&url);

		save(config.get(), destStorage, format, handler);

		includeStack.pop_back();

		setName(format.getConfigName(destStorage));
		setUrl(destStorage.getUrl());
	}
	RETHROW_AS_IO_ERROR(IoError::Save, destStorage.getUrl());
}

/**
 * Saves the configuration document in its original storage(s).
 * The information about where to save each subtree comes from the previous
 * call to any "load" function, which initializes the "includeURL" property of
 * configuration objects.
 * Nevertheless, it can also be modified using
 * CCfg::Object::setIncludeUrl(const std::string&), or created from scratch
 * using first  CCfg::Io::Document::setConfig(Object*) to associate a
 * configuration to the document, and then
 * CCfg::Object::setIncludeUrl(const std::string&).
 * In this case, it is worth to note that the configuration root must have an
 * includeURL otherwise, nothing will be saved.
 * @throws Document::IoError reports that the document cannot be saved.
 */
void Document::save()
{
	assert(config.get() != 0); // TODO: throw exception instead
	try
	{
		Url url(getUrl());
		includeStack.clear();
		includeStack.push_back(&url);

		saveToUrl(config.get(), url, *this);

		includeStack.pop_back();
	}
	RETHROW_AS_IO_ERROR(IoError::Save, getUrl());
}

/**
 * Includes a configuration into a destination object.
 * Relative URLs must be relative to the include directory of the destination
 * object.
 */
void Document::include(Object& destination, const Url& resourceUrl)
{
	Url parentIncludeUrl = computeQualifiedIncludeUrl(*destination.getIncludeRoot());
	Url url = resourceUrl;
	toContextRelativeUrl(url, parentIncludeUrl);

	std::auto_ptr<Object> includeRoot = loadFromUrl(url, *this);
	destination.import(includeRoot);
	destination.setIncludeUrl(resourceUrl.format());
}

/**
 * Imports a configuration into a destination object.
 * Relative URLs must be relative to the include directory of the destination
 * object.
 */
void Document::import(Object& destination, const Url& resourceUrl)
{
	Url parentIncludeUrl = computeQualifiedIncludeUrl(*destination.getIncludeRoot());
	Url url = resourceUrl;
	toContextRelativeUrl(url, parentIncludeUrl);

	std::auto_ptr<Object> importedRoot = loadFromUrl(url, *this);
	destination.import(importedRoot);
}

/**
 * Exports the an Object and its subtree to an URL.
 */
void Document::exportSubConfig(const Object& root, const Url& resourceUrl) const
{
	saveToUrl(&root, resourceUrl, NoOpIncludeHandler::getInstance());
}

/**
 * Resolves the include directive in an object.
 * The includeURL attribute must be defined.
 * @note If the object is not empty it will be cleared.
 * Includes will be processed recursively.
 * @param includerObject the object with an include directive.
 * @throws Storage::OpenError reports there was an error opening the include
 * storage.
 * @throws Format::SyntaxError reports the document cannot be parsed because of a
 * syntax error.
 * @throws Format::SemanticError reports the included document cannot be parsed
 * because the semantics constraints on the format are not respected.
 * @throws Factory::UnhandledUrl reports we don't know how to open the resource
 * designated by an include Url.
 * @throws Document::RecursiveInclude reports an infinite recursive inclusion cycle.
 */
void Document::parseInclude(Object& includerObject)
{
	includerObject.clear(); // ensure the object is empty

	Url url(includerObject.getIncludeUrl());
	toContextRelativeUrl(url, *includeStack.back());

	if(isRecursiveInclude(url))
	{
		throw RecursiveInclude(url.format(), includeStack.back()->format());
	}

	includeStack.push_back(&url);
	std::auto_ptr<Object> includeRoot = loadFromUrl(url, *this);
	includeStack.pop_back();

	includerObject.import(includeRoot);
}

/**
 * Saves an included subtree in the storage designated by its Url.
 * @param includerObject the object with an include directive.
 * @throws Storage::OpenError reports there was an error opening the storage.
 */
bool Document::formatInclude(const Object& includerObject)
{
	Url url(includerObject.getIncludeUrl());
	toContextRelativeUrl(url, *includeStack.back());

	includeStack.push_back(&url);
	saveToUrl(&includerObject, url, *this);
	includeStack.pop_back();
	return true;
}

/**
 * Loads a configuration from a storage without changing document state.
 * @param srcUrl the URL of the storage.
 * @param handler the callback object that will have to process includes.
 * @throws Storage::OpenError reports there was an error opening the storage.
 * @throws Format::SyntaxError reports the document cannot be parsed because of
 * a syntax error.
 * @throws Format::SemanticError reports the document cannot be parsed because the
 * semantics constraints on the format are not respected.
 * @throws Factory::UnhandledUrl reports we don't know how to open the resource
 * designated by the Url.
 * @throws Document::RecursiveInclude reports an infinite recursive inclusion cycle.
 */
std::auto_ptr<Object> Document::loadFromUrl(const Url& srcUrl,
	IncludeHandler& handler)
{
	std::auto_ptr<Format> format = FormatFactory::instance().createFormat(srcUrl);
	std::auto_ptr<Storage> storage = StorageFactory::instance().createStorage(srcUrl);

	return load(*storage, *format, handler);
}

/**
 * Saves a configuration in a storage.
 * @param config the configuration to save.
 * @param destUrl the URL of the storage.
 * @param handler the callback object that will have to process includes.
 * @throws Storage::OpenError reports there was an error opening the storage.
 * @throws Factory::UnhandledUrl reports we don't know how to open the resource
 * designated by the Url.
 */
void Document::saveToUrl(const Object *config,
	const Url& destUrl,
	IncludeHandler& handler) const
{
	Url url(destUrl);
	std::auto_ptr<Format> format = FormatFactory::instance().createFormat(destUrl);
	std::auto_ptr<Storage> storage = StorageFactory::instance().createStorage(destUrl);

	save(config, *storage, *format, handler);
}

/**
 * Loads a configuration from a storage without changing document state.
 * The root object automatically gets the URL of the storage where it comes
 * from.
 * @param srcStorage the storage where the configuration should be read.
 * @param format the configuration decoding format.
 * @param handler the callback object that will have to process includes.
 * @throws Storage::OpenError reports there was an error opening the storage.
 * @throws Format::SyntaxError reports the document cannot be parsed because of a
 * syntax error.
 * @throws Format::SemanticError reports the document cannot be parsed because the
 * semantics constraints on the format are not respected.
 * @throws Url::ParseError reports an include URL was invalid.
 * @throws Factory::UnhandledUrl reports we don't know how to open the resource
 * designated by an include Url.
 * @throws Document::RecursiveInclude reports an infinite recursive inclusion cycle.
 */
std::auto_ptr<Object> Document::load(Storage& srcStorage,
	Format& format,
	IncludeHandler& handler)
{
	return format.parse(srcStorage, handler);
}

/**
 * Saves a configuration in a storage.
 * @param config the configuration to save.
 * @param destStorage the destination storage.
 * @param format the configuration encoding format.
 * @param handler the callback object that will have to process includes.
 * @throws Storage::OpenError reports there was an error opening the storage.
 * @throws Url::ParseError reports an include URL was invalid.
 * @throws Factory::UnhandledUrl reports we don't know how to open the resource
 * designated by an include Url.
 */
void Document::save(const Object *config,
	Storage& destStorage,
	Format& format,
	IncludeHandler& handler) const
{
	format.format(config, destStorage, handler);
}

/**
 * Check whether the specified include Url produces a recursive include.
 * @param url the Url to check.
 * @return Returns true if the include is recursive, false otherwise.
 */
bool Document::isRecursiveInclude(const Url& url) const
{
	for(UrlRefStack::const_iterator iUrl = includeStack.begin();
		iUrl != includeStack.end(); ++iUrl)
	{
		if(url.isSameResource(**iUrl))
		{
			return true;
		}
	}
	return false;
}

/**
 * Getter of the "config" property.
 * @return Returns the configuration resulting of a preceding call to a public
 * load method.
 */
const Object* Document::getConfig() const
{
	return config.get();
}

/**
 * Getter of the "config" property.
 * @return Returns the configuration tied to this document.
 * It results of a preceding call to a public load or setConfig method.
 */
Object* Document::getConfig()
{
	return config.get();
}

/**
 * Setter of the "config" property.
 * @param newCfg the new configuration tied to this document.
 */
void Document::setConfig(std::auto_ptr<Object> newCfg)
{
	config = newCfg;
}

/**
 * Release ownership of the config tree associated to this document.
 * Only load functions or setConfig can be called afterwards.
 */
std::auto_ptr<Object> Document::detachConfig()
{
	return config;
}

/**
 * Returns the URL where the document root will be saved.
 */
const std::string& Document::getUrl() const
{
	return config->getIncludeUrl();
}

/**
 * Defines the URL where the document root will be saved.
 */
void Document::setUrl(const std::string& url)
{
	config->setIncludeUrl(url);
}

/**
 * Returns the name of this document.
 */
const std::string& Document::getName() const
{
	return name;
}

/**
 * Defines the name of this document.
 */
void Document::setName(const std::string& newName)
{
	name = newName;
}

/**
 * Returns the qualified include URL of the specified attribute.
 * By qualified URL we mean: relative to the current directory, or absolute.
 * If there are several level of includes using relative path, computing the
 * include path in the inner level requires some computations.
 * @param obj the object of which to retrieve the include path.
 * @return Returns an empty URL if the include URL of obj is empty, returns the
 * qualified URL otherwise.
 */
Url Document::computeQualifiedIncludeUrl(const Object& obj)
{
	Url qualifiedUrl;
	const std::string& includeUrl = obj.getIncludeUrl();
	if(not includeUrl.empty())
	{
		qualifiedUrl.parse(includeUrl);

		const Attribute* current = &obj;
		const Attribute* includeRoot = obj.getIncludeRoot();

		while(current != includeRoot and qualifiedUrl.isRelative())
		{
			current = includeRoot;
			if(qualifiedUrl.isRelative())
			{
				Url url(current->getIncludeUrl());
				toContextRelativeUrl(qualifiedUrl, url);
			}

			includeRoot = current->getIncludeRoot();
		}
	}

	return qualifiedUrl;
}

/**
 * Converts an URL to a context relative URL in the include hierarchy, if necessary.
 * @param url the URL that will be transformed to a context relative URL.
 * @param parentInclude the URL of the parent configuration where the first URL
 * is defined, and to which it is relative.
 */
Url& Document::toContextRelativeUrl(Url& url, const Url& parentInclude)
{
	if(url.isRelative())
	{
		std::string qualifiedPath = parentInclude.getDirectory();
		if(not qualifiedPath.empty()) // if the parent is not a file in the current directory
		{
			qualifiedPath += Url::PATH_SEPARATOR;
		}
		qualifiedPath += url.getPath();
		url.setPath(qualifiedPath);
	}
	return url;
}



///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Document::IoError class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param ioType the type of IO operation that failed.
 * @param docUrl the URL of the document for which IO operation failed.
 * @param cause the message describing the cause of the error.
 */
Document::IoError::IoError(Document::IoError::IoType ioType,
	const std::string& docUrl,
	const std::string& cause)
	: Exception(), ioType(ioType), docUrl(docUrl), cause(cause)
{
}

/**
 * Destructor.
 */
Document::IoError::~IoError() throw()
{
}

/**
 * Getter of the "ioType" property.
 * @return Returns the type of IO operation that failed.
 */
Document::IoError::IoType Document::IoError::getIoType() const
{
	return ioType;
}

/**
 * Getter of the "docUrl" property.
 * @return Returns the URL of the document for which IO operation failed.
 */
const std::string& Document::IoError::getDocUrl() const
{
	return docUrl;
}

/**
 * Getter of the "cause" property.
 * @return Returns the message describing the cause of the error.
 */
const std::string& Document::IoError::getCause() const
{
	return cause;
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string Document::IoError::formatMessage() const
{
	std::ostringstream message;
	message << "Error " << (ioType == Load ? "loading" : "saving") << " the document";
	if(not docUrl.empty())
	{
		message << " located at: \"" <<  docUrl << '\"';
	}
	message << '.' << std::endl << cause;
	return message.str();
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Document::RecursiveInclude class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param recursiveUrl the URL that already exists in the include stack and
 * that would generate an infinite recursion.
 * @param parentFileUrl the URL of the file where the bad include instruction
 * is.
 */
Document::RecursiveInclude::RecursiveInclude(const std::string& recursiveUrl,
	const std::string& parentFileUrl)
	: recursiveUrl(recursiveUrl), parentFileUrl(parentFileUrl)
{
}

/**
 * Destructor.
 */
Document::RecursiveInclude::~RecursiveInclude() throw()
{
}

/**
 * Getter of the "recursiveUrl" property.
 * @return Returns the URL that would generate an infinite recursion.
 */
const std::string& Document::RecursiveInclude::getRecursiveUrl() const
{
	return recursiveUrl;
}

/**
 * Getter of the "parentFileUrl" property.
 * @return Returns the URL of the file where the bad include instruction is.
 */
const std::string& Document::RecursiveInclude::getParentFileUrl() const
{
	return parentFileUrl;
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string Document::RecursiveInclude::formatMessage() const
{
	std::ostringstream message;
	message << "The file \"" << parentFileUrl <<
		"\" cannot include the file \"" << recursiveUrl <<
		"\" because this would begin an infinite recursion.";
	return message.str();
}

}

}
