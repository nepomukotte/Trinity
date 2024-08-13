// =============================================================================
// Url.h
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

#ifndef CCfg_Io_Url_h_INCLUDED
#define CCfg_Io_Url_h_INCLUDED

#include "../Exception.h"
#include <string>
#include <map>
#include <ostream>

namespace CCfg
{

namespace Io
{
class UrlParser;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class enable URL (Uniform Resource Locator) parsing and processing.
 * The parser only implements a subset of the specifications:
 * - generic URLs: scheme://user:password\@host:port/path?query.<br>
 *   file urls fall into this category without distinction.
 * - relative URLs segment/segment/segment etc...
 * - other URLs (mail, newsgroups, etc...) are not supported.
 *
 * @author fchateau
 */
class Url
{
protected:
	typedef std::map<std::string, std::string> QueryMap;

	/// the scheme, e.g.: file, http, ftp, etc. There is not restriction on the
	/// scheme as long as it respect the grammar.
	std::string scheme;
	std::string user;	///< the user name for connecting the host.
	std::string password;	///< the password for connecting the host.
	std::string host;	///< the host dns (domains are supported), or ipv4 address.
	int port;			///< the port to connect on.
	std::string path;	///< the resource path on the host, or the relative path (if relative URL).
	QueryMap query;		///< the URL query made of pair key=value beginning by '?' and separated by '&'

public:
	friend class UrlParser;
	class ParseError;
	struct PathNormalizer;

	static const char PATH_SEPARATOR = '/';
	static const char FILE_SCHEME[];

	// construction/destruction
	Url();
	explicit Url(const std::string& url);
	Url(const Url& parser);
	~Url();
	Url& operator=(const Url& parser);

	// string conversions.
	void parse(const std::string& url);
	std::string format() const;

	// getters
	const std::string& getScheme() const;
	const std::string& getUser() const;
	const std::string& getPassword() const;
	const std::string& getHost() const;
	int getPort() const;
	const std::string& getPath() const;
	const QueryMap& getQuery() const;
	std::string getQueryItem(const std::string& key) const;

	// helper getters
	bool isEmpty() const;
	bool isFile() const;
	bool isRelative() const;
	bool hasScheme() const;
	bool hasUser() const;
	bool hasPassword() const;
	bool hasHost() const;
	bool hasPath() const;
	bool hasQuery() const;
	bool hasQueryItem(const std::string& key);
	std::string getDirectory() const;
	std::string getFile() const;
	std::string getExtension() const;

	// setters
	void setScheme(const std::string& value);
	void setUser(const std::string& value);
	void setPassword(const std::string& value);
	void setHost(const std::string& value);
	void setPort(int value);
	void setPath(const std::string& value);
	void setQuery(const QueryMap&  value);

	// url alteration
	void addQueryItem(const std::string& key, const std::string& value);
	void removeQueryItem(const std::string& key);
	void appendPath(const std::string& pathElement);

	// comparison
	bool operator==(const Url& rvalue);
	bool operator!=(const Url& rvalue);
	bool isSameResource(const Url& url) const;

	void reset();

	static Url fromPath(const std::string& path);
	std::string toPath() const;

	friend std::ostream& operator<<(std::ostream& stream, const Url& url);
};

std::ostream& operator<<(std::ostream& stream, const Url& url);


/**
 * This class is a stream manipulator that handle conversion between file paths
 * and URL paths.
 * In file paths most characters are allowed like spaces, angled bracked...,
 * whereas URL paths have more restrictions.
 *
 * @author fchateau
 */
struct Url::PathNormalizer
{
	static const char ALLOWED_SPECIAL_CHARACTERS[];

	std::string& path;
	PathNormalizer(std::string& path);
	PathNormalizer(const PathNormalizer& rvalue);
	~PathNormalizer();

	static bool isAllowedSpecial(char c);
	static bool isAllowedCharacter(char c);

	static std::string escape(const std::string& str);
	static std::string unescape(const std::string& str);

private:
	PathNormalizer();
	PathNormalizer& operator=(const PathNormalizer& rvalue);
};

std::ostream& operator<<(std::ostream& stream, const Url::PathNormalizer& normalizer);
std::istream& operator>>(std::istream& stream, const Url::PathNormalizer& normalizer);

/**
 * This class reports an error while parsing a URL.
 *
 * @author fchateau
 */
class Url::ParseError : public Exception
{
	std::string url;	///< the erroneous URL.
	int errorCursor;	///< the position of the character where parsing failed.

public:
	// construction/destruction
	ParseError(const std::string& url, int errorCursor);
	virtual ~ParseError() throw();

	// getters
	const std::string& getUrl() const;
	int getErrorCursor() const;

private:
	std::string formatMessage() const;
};





///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Url inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Getter of the "scheme" property.
 * @return Returns the URL scheme.
 */
inline const std::string& Url::getScheme() const
{
	return scheme;
}

/**
 * Getter of the "user" property.
 * @return Returns the URL user.
 */
inline const std::string& Url::getUser() const
{
	return user;
}

/**
 * Getter of the "password" property.
 * @return Returns the URL password.
 */
inline const std::string& Url::getPassword() const
{
	return password;
}

/**
 * Getter of the "host" property.
 * @return Returns the URL host.
 */
inline const std::string& Url::getHost() const
{
	return host;
}

/**
 * Getter of the "port" property.
 * @return Returns the URL port.
 */
inline int Url::getPort() const
{
	return port;
}

/**
 * Getter of the "path" property.
 * @return Returns the URL path.
 */
inline const std::string& Url::getPath() const
{
	return path;
}

/**
 * Getter of the "query" property.
 * @return Returns the URL query.
 */
inline const Url::QueryMap& Url::getQuery() const
{
	return query;
}

/**
 * Returns the value of the specified query item.
 */
inline std::string Url::getQueryItem(const std::string& key) const
{
	std::string value;
	QueryMap::const_iterator item = query.find(key);
	if(item != query.end())
	{
		value = item->second;
	}
	return value;
}

/**
 * Determines whether the URL is empty (all relevant fields are empty).
 */
inline bool Url::isEmpty() const
{
	return scheme.empty() and host.empty() and path.empty();
}

/**
 * Determines whether the URL is relative or absolute.
 * @return Return true if the URL is relative, false otherwise.
 */
inline bool Url::isRelative() const
{
	return scheme.empty() and (path.empty() or path[0] != Url::PATH_SEPARATOR);
}

/**
 * Determines whether the URL represents a file resource.
 */
inline bool Url::isFile() const
{
	return (scheme.empty() or scheme == FILE_SCHEME) and not path.empty() and host.empty();
}

/**
 * Determines whether the URL contains a scheme.
 */
inline bool Url::hasScheme() const
{
	return not scheme.empty();
}

/**
 * Determines whether the URL contains a user.
 * @return Returns true if the URL has a user.
 */
inline bool Url::hasUser() const
{
	return not user.empty();
}

/**
 * Determines whether the URL contains a password.
 * @return Returns true if the URL has a password.
 */
inline bool Url::hasPassword() const
{
	return not password.empty();
}

/**
 * Determines whether the URL contains a host.
 * @return true if the URL has a host.
 */
inline bool Url::hasHost() const
{
	return not host.empty();
}

/**
 * Determines whether the URL contains a path.
 * @return Returns true if the URL has a path.
 */
inline bool Url::hasPath() const
{
	return not path.empty();
}

/**
 * Determines whether the URL contains a query.
 * @return Returns true if the URL has a query.
 */
inline bool Url::hasQuery() const
{
	return not query.empty();
}

/**
 * Determines whether the URL query has a value for the specified key.
 * @param key the key whose existence is to test.
 * @return Returns true if the key exists in the query.
 */
inline bool Url::hasQueryItem(const std::string& key)
{
	return query.find(key) != query.end();
}

/**
 * Setter of the "scheme" property.
 * @param value the new scheme.
 */
inline void Url::setScheme(const std::string& value)
{
	scheme = value;
}

/**
 * Setter of the "user" property.
 * @param value the new user.
 */
inline void Url::setUser(const std::string& value)
{
	user = value;
}

/**
 * Setter of the "password" property.
 * @param value the new password.
 */
inline void Url::setPassword(const std::string& value)
{
	password = value;
}

/**
 * Setter of the "host" property.
 * @param value the new host.
 */
inline void Url::setHost(const std::string& value)
{
	host = value;
}

/**
 * Setter of the "port" property.
 * @param value the new port.
 */
inline void Url::setPort(int value)
{
	port = value;
}

/**
 * Setter of the "path" property.
 * @param value the new path.
 */
inline void Url::setPath(const std::string& value)
{
	path = value;
}

/**
 * Setter of the "scheme" property.
 * @param value the new scheme.
 */
inline void Url::setQuery(const QueryMap&  value)
{
	query = value;
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Url::PathNormalizer inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
inline Url::PathNormalizer::PathNormalizer(std::string& path)
	: path(path)
{
}

/**
 * Copy constructor.
 */
inline Url::PathNormalizer::PathNormalizer(const PathNormalizer& rvalue)
	: path(rvalue.path)
{
}

/**
 * Destructor.
 */
inline Url::PathNormalizer::~PathNormalizer()
{
}

}
}

#endif // CCfg_Io_Url_h_INCLUDED
