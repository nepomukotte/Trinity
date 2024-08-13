// =============================================================================
// Document.h
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

#ifndef CCfg_Io_Document_h_INCLUDED
#define CCfg_Io_Document_h_INCLUDED

#include "IncludeHandler.h"
#include "../Exception.h"
#include <string>
#include <memory>
#include <vector>

namespace CCfg
{
// forward declarations
class Object;
class Attribute;

namespace Io
{
// forward declarations
class Format;
class Storage;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is responsible for loading and saving configurations in any format.
 * It owns the configuration root object, so its destruction will free the entire
 * configuration tree.
 * New format and storage may be defined by implementing CCfg::Io::Format,
 * CCfg::Io::Storage, and providing a factory implementing CCfg::Io::Factory
 * to instanciate these classes according to an Url.
 * The document factory must then be redefined to register the extensions. The
 * default factory is defined in the constructor, subclassing is necessary to
 * change the default factory.
 *
 * @author fchateau
 */
class Document : public IncludeHandler
{
	std::auto_ptr<Object> config;	///< the configuration tied to this document.
	std::string name;	///< the name of the document.

	typedef std::vector<Url*> UrlRefStack;

	// transient attributes
	/// the stack of URLs of nested include files while loading or saving.
	/// Warning: pointers reference automatic objects created on the
	/// application stack through the recursive process.
	UrlRefStack includeStack;

public:
	// inner classes declaration
	class IoError;
	class RecursiveInclude;

	// @{ construction/destruction
	Document();
	explicit Document(const std::string& path);
	explicit Document(const Url& srcUrl);
	virtual ~Document();
	// @}

	// @{ document loading and saving
	Object& loadFromFile(const std::string& path, bool ignoreIncludes = false);
	void saveToFile(const std::string& path, bool flatten = true);
	Object& loadFromUrl(const Url& srcUrl, bool ignoreIncludes = false);
	void saveToUrl(const Url& destUrl, bool flatten = true);
	Object& load(Storage& srcStorage, Format& format, bool ignoreIncludes = true);
	void save(Storage& destStorage, Format& format, bool flatten = true);
	void save();
	// @}

	// @{ document import/export
	void include(Object& destination, const Url& resourceUrl);
	void import(Object& destination, const Url& resourceUrl);
	void exportSubConfig(const Object& root, const Url& resourceUrl) const;
	// @}

	// @{ getters and setters
	const Object* getConfig() const;
	Object* getConfig();
	void setConfig(std::auto_ptr<Object> newCfg);
	std::auto_ptr<Object> detachConfig();

	const std::string& getUrl() const;
	const std::string& getName() const;
	void setUrl(const std::string& url);
	void setName(const std::string& newName);
	// @}

	static Url computeQualifiedIncludeUrl(const Object& obj);
	static Url& toContextRelativeUrl(Url& url, const Url& parentInclude);

protected:
	// include callbacks
	void parseInclude(Object& includerObject);
	bool formatInclude(const Object& includerObject);

	// internal loading and saving method
	std::auto_ptr<Object> loadFromUrl(const Url& srcUrl, IncludeHandler& handler);
	void saveToUrl(const Object *config, const Url& destUrl, IncludeHandler& handler) const;

	std::auto_ptr<Object> load(Storage& srcStorage, Format& format, IncludeHandler& handler);
	void save(const Object *config, Storage& destStorage, Format& format, IncludeHandler& handler) const;

	bool isRecursiveInclude(const Url& url) const;

private:
	// non copyable semantics
	Document(const Document&);
	Document& operator=(const Document&);
};


/**
 * This exception reports that an error occured while loading or saving a
 * document.
 */
class Document::IoError : public Exception
{
public:
	enum IoType {Load, Save};
private:
	IoType ioType;		///< the type of IO operation that failed.
	std::string docUrl;	///< the URL of the document for which IO operation failed.
	std::string cause;	///< the message describing the cause of the error.

public:
	// construction/destruction
	IoError(IoType ioType, const std::string& docUrl, const std::string& cause);
	~IoError() throw();

	IoType getIoType() const;
	const std::string& getDocUrl() const;
	const std::string& getCause() const;

private:
	std::string formatMessage() const;
};


/**
 * This exception reports that there is direct or indirect recursivity in
 * include directives.
 */
class Document::RecursiveInclude : public Exception
{
	std::string recursiveUrl;	///< the URL that would generate an infinite recursion.
	std::string parentFileUrl;	///< the URL of the file where the bad include instruction is.

public:
	RecursiveInclude(const std::string& recursiveUrl, const std::string& parentFileUrl);
	~RecursiveInclude() throw();

	const std::string& getRecursiveUrl() const;
	const std::string& getParentFileUrl() const;

private:
	std::string formatMessage() const;
};

}

}

#endif // CCfg_Io_Document_h_INCLUDED
