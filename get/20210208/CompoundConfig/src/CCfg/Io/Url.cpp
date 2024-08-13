// =============================================================================
// Url.cpp
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

#include "Url.h"
#include <boost/version.hpp>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/bind.hpp>
#include <bitset>
#include <sstream>
#include <iomanip>
#include <cstdlib>

namespace boost_spirit_classic = boost::spirit::classic;
using std::string;
using namespace boost;
using namespace boost::spirit::classic;

namespace CCfg
{

namespace Io
{
///////////////////////////////////////////////////////////////////////////////
// private classes declarations
///////////////////////////////////////////////////////////////////////////////

class UrlParser
{
	Url* url;

	enum UrlValidations {SCHEME, AUTH, HOSTPORT, LOGIN, PATH, QUERY};
	std::bitset<8> valid;

public:
	UrlParser(Url *url);

	void onScheme(const char* begin, const char* end);
	void onUser(const char* begin, const char* end);
	void onPassword(const char* begin, const char* end);
	void onHost(const char* begin, const char* end);
	void onPort(int value);
	void onPath(const char* begin, const char* end);
	void onQueryItem(const char* begin, const char* end);

	void validateScheme(const char*, const char*);
	void validateAuth(const char*, const char*);
	void validateHostport(const char*, const char*);
	void validateLogin(const char*, const char*);
	void validatePath(const char*, const char*);
	void validateQuery(const char*, const char*);
	void validateUrl(const char*, const char*);
};

struct UrlGrammar : grammar<UrlGrammar>
{
	UrlGrammar(UrlParser* parser) : parser(parser)
	{
	}

	template <typename ScannerT>
	struct definition
	{
		definition(const UrlGrammar& self)
		{
			#define BIND(fun)	bind(&UrlParser::fun, self.parser, _1, _2)
			#define BIND1(fun)	bind(&UrlParser::fun, self.parser, _1)

			url =
			(
				url_main	= 	(absolute_url | relative_url)[BIND(validateUrl)],
				absolute_url= 	scheme_spec[BIND(validateScheme)] >>
								!(login[BIND(validateLogin)]) >>
								!(path_spec[BIND(validatePath)]) >>
								!(query_spec[BIND(validateQuery)]),
				relative_url=	path[BIND(onPath)][BIND(validatePath)],
				scheme_spec	=	scheme[BIND(onScheme)] >> "://",
				path_spec	=	(Url::PATH_SEPARATOR >> path[BIND(onPath)]),
				query_spec	=	('?' >> query),
				scheme		=	+(lower_p | digit_p | '+' | '-' | '.'),
				login		=	!auth[BIND(validateAuth)] >>
								hostport[BIND(validateHostport)]
			);
			query =
			(
				query_main		=	fieldspec[BIND(onQueryItem)] >>
									*('&' >> fieldspec[BIND(onQueryItem)]),
				fieldspec		=	fieldname >> '=' >> fieldvalue,
				fieldsymbol		=	*(uchar | ':' | ';' | '@')
				// fieldname	=	fieldsymbol
				// fieldvalue	= 	fieldsymbol
			);
			path =
			(
				path_main	=	segment >> *(Url::PATH_SEPARATOR >> segment),
				segment		=	*(uchar | ':' | '@' | '&' | '=')
			);
			hostport =
			(
				hostport_main	=	host[BIND(onHost)] >> !(':' >> int_p[BIND1(onPort)]),
				host			=	hostname | hostipv4,
				hostname 		=	*(domainlabel >> '.') >> toplabel,
				domainlabel		=	alnum_p >> *(*ch_p('-') >> alnum_p),
				toplabel		=	alpha_p >> *(*ch_p('-') >> alnum_p),
				hostipv4		=	int_p >> '.' >> int_p >> '.' >> int_p >> '.' >> int_p
			);
			auth =
			(
				auth_main	=	user[BIND(onUser)] >>
								!(':' >> password[BIND(onPassword)]) >>
								'@',
				authsymbol	=	*(uchar | ';' | '?' | '&' | '=')
				// user		=	authsymbol
				// password	=	authsymbol
			);
			uchar =
			(
				uchar_main	=	unreserved | escape,
				unreserved	=	alpha_p | digit_p | safe | extra,
				safe		=	ch_p('$') | '-' | '_' | '.' | '+',
				extra		=	ch_p('!') | '*' | '\'' | '(' | ')' | ',',
				escape		=	ch_p('%') >> xdigit_p >> xdigit_p
			);
		}

		const rule<ScannerT>& start() const
		{
			return url;
		}

		rule<ScannerT> url;
		subrule<0> url_main;
		subrule<1> absolute_url;
		subrule<2> relative_url;
		subrule<3> scheme_spec;
		subrule<4> path_spec;
		subrule<5> query_spec;
		subrule<6> scheme;
		subrule<7> login;


		rule<ScannerT> query;
		subrule<0> query_main;
		subrule<1> fieldspec;
		subrule<2> fieldsymbol;
		subrule<2> fieldname;	// alias for fieldsymbol
		subrule<2> fieldvalue;	// alias for fieldsymbol

		rule<ScannerT> path;
		subrule<0> path_main;
		subrule<1> segment;

		rule<ScannerT> hostport;
		subrule<0> hostport_main;
		subrule<1> host;
		subrule<2> hostname;
		subrule<3> domainlabel;
		subrule<4> toplabel;
		subrule<5> hostipv4;

		rule<ScannerT> auth;
		subrule<0> auth_main;
		subrule<1> authsymbol;
		subrule<1> user;		// alias for authsymbol
		subrule<1> password;	// alias for authsymbol

		rule<ScannerT> uchar;
		subrule<0> uchar_main;
		subrule<1> unreserved;
		subrule<2> safe;
		subrule<3> extra;
		subrule<4> escape;
	};

	UrlParser *parser;
};





///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Url class methods body
///////////////////////////////////////////////////////////////////////////////

const char Url::PATH_SEPARATOR;
const char Url::FILE_SCHEME[] = "file";

/**
 * Default constructor.
 */
Url::Url()
	: scheme(),
	user(),
	password(),
	host(),
	port(0),
	path(),
	query()
{
}

/**
 * Transtyping constructor. Converts a string to a URL by parsing it.
 */
Url::Url(const std::string& url)
{
	parse(url);
}

/**
 * Copy constructor. Make a copy of an URL.
 * @param url the URL to copy.
 */
Url::Url(const Url& url)
	: scheme(url.scheme),
	user(url.user),
	password(url.password),
	host(url.host),
	port(url.port),
	path(url.path),
	query(url.query)
{
}

/**
 * Destructor.
 */
Url::~Url()
{
}

/**
 * Assignment operator.
 * @param url the URL to copy.
 */
Url& Url::operator=(const Url& url)
{
	scheme = url.scheme;
	user = url.user;
	password = url.password;
	host = url.host;
	port = url.port;
	path = url.path;
	query = url.query;
	return *this;
}

/**
 * Initializes an URL from a file path.
 * Special characters like spaces are escaped on the process to make the path
 * compliant with URL specifications.
 */
Url Url::fromPath(const std::string& path)
{
	Url url;
	url.setPath(PathNormalizer::escape(path));
	return url;
}

/**
 * Converts this URL to a file path.
 * Escaped characters are restored.
 */
std::string Url::toPath() const
{
	return PathNormalizer::unescape(path);
}

/**
 * Resets *this with the parsing of a string representing a URL.
 * @param url the string URL.
 * @throws Url::ParseError reports that the URL was invalid.
 */
void Url::parse(const std::string& url)
{
	reset();
	UrlParser parser(this);
	UrlGrammar g(&parser);
	const parse_info<const char*>& ret = boost_spirit_classic::parse(url.c_str(), g);
	if(not ret.full or not ret.hit)
	{
		throw ParseError(url, ret.length);
	}
}

/**
 * Convert the URL to its string representation.
 * @return Returns the string representation of the URL.
 */
std::string Url::format() const
{
	std::ostringstream url;
	url << *this;
	return url.str();
}


/**
 * Helper to get the directory part of the path.
 * @return Returns the directory part of the path (before the last slash).
 */
std::string Url::getDirectory() const
{
	size_t idx = path.find_last_of(Url::PATH_SEPARATOR);
	return (idx == std::string::npos) ? "" : path.substr(0, idx);
}

/**
 * Helper to get the file part of the path.
 * @return Returns the file part of the path (after the last slash).
 */
std::string Url::getFile() const
{
	return path.substr(path.find_last_of(Url::PATH_SEPARATOR) + 1);
	// works also if there is no slash (because it return -1, so we get the
	// whole string).
}

/**
 * Helper to get the extension part of the path.
 * @return Returns the extension part of the path (after the last dot).
 */
std::string Url::getExtension() const
{
	size_t pos = path.find_last_of('.');
	if(pos != std::string::npos)
	{
		return path.substr(pos + 1);
	}
	else
	{
		return std::string();
	}
}

/**
 * Adds a new item (couple (key, value)) in the URL query.
 * @param key the item key.
 * @param value the item value.
 */
void Url::addQueryItem(const std::string& key, const std::string& value)
{
	query.insert(make_pair(key, value));
}

/**
 * Removes a query item.
 * @param key the key of the item to remove.
 */
void Url::removeQueryItem(const std::string& key)
{
	query.erase(key);
}

/**
 * Appends the specified path element to the path.
 * @note a path separator is automatically appended before if the path is not empty.
 */
void Url::appendPath(const std::string& pathElement)
{
	if(not path.empty())
	{
		path += PATH_SEPARATOR;
	}
	path += pathElement;
}

/**
 * Equality comparison operator.
 */
bool Url::operator==(const Url& rvalue)
{
	return scheme == rvalue.scheme and
		user == rvalue.user and
		password == rvalue.password and
		host == rvalue.host and
		port == rvalue.port and
		path == rvalue.path and
		query == rvalue.query;
}

/**
 * Unequality comparison operator.
 */
bool Url::operator!=(const Url& rvalue)
{
	return !operator==(rvalue);
}

/**
 * Tells whether two URLs undoublty reference the same resource.
 * @note
 * This function only compares URL fields, it does not try to see if different URLs
 * reference the same physical resource. The authentication informations are ignored
 * @param url the URL to compare to *this.
 * @return Returns true if URLs have exactly identical location informations,
 * false otherwise.
 */
bool Url::isSameResource(const Url& url) const
{
	return scheme == url.getScheme() and
		host == url.getHost() and
		port == url.getPort() and
		path == url.getPath() and
		query == url.getQuery();
}

/**
 * Resets the URL (every property becomes undefined).
 */
void Url::reset()
{
	scheme.clear();
	user.clear();
	password.clear();
	host.clear();
	port = 0;
	path.clear();
	query.clear();
}

/**
 * Stream insertion operator. Outputs the URL on a stream.
 * @param stream the stream where the URL will be appended.
 * @param url the URL to write.
 * @return Returns the altered stream.
 */
std::ostream& operator<<(std::ostream& stream, const Url& url)
{
	if(not url.hasScheme())
	{
		stream << url.getPath();
	}
	else
	{
		stream << url.getScheme() << "://";
		if(url.hasUser())
		{
			stream << url.getUser();
			if(url.hasPassword())
			{
				stream << ':' << url.getPassword();
			}
			stream  << '@';
		}
		if(url.hasHost())
		{
			stream << url.getHost();
			if(url.getPort() != 0)
			{
				stream << ':' << url.getPort();
			}
		}
		if(url.hasPath())
		{
			stream << Url::PATH_SEPARATOR << url.getPath();
		}
		if(url.hasQuery())
		{
			stream << '?';
			const Url::QueryMap& query = url.getQuery();
			Url::QueryMap::const_iterator entry = query.begin();
			while(entry != query.end())
			{
				stream << entry->first << '=' << entry->second;
				if( ++entry != query.end())
				{
					stream << '&';
				}
			}
		}
	}
	return stream;
}


////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Url::PathNormalizer class methods body
////////////////////////////////////////////////////////////////////////////////

const char Url::PathNormalizer::ALLOWED_SPECIAL_CHARACTERS[] = "/$-_.+!*'(),:@&=";

/**
 * Converts a file path to an URL path.
 */
std::string Url::PathNormalizer::escape(const std::string& str)
{
	std::ostringstream escapedStr;
	escapedStr << PathNormalizer(const_cast<std::string&>(str));
	return escapedStr.str();
}

/**
 * Converts an URL path to a file path.
 */
std::string Url::PathNormalizer::unescape(const std::string& str)
{
	std::string unescapedStr;
	std::istringstream escapedStr(const_cast<std::string&>(str));
	escapedStr >> PathNormalizer(unescapedStr);
	return unescapedStr;
}

/**
 * Returns whether the specified character is part of the authorized special
 * characters in URL paths.
 */
bool Url::PathNormalizer::isAllowedSpecial(char c)
{
	return
		std::find(&ALLOWED_SPECIAL_CHARACTERS[0],
			&ALLOWED_SPECIAL_CHARACTERS[sizeof ALLOWED_SPECIAL_CHARACTERS],
			c)
		!= &ALLOWED_SPECIAL_CHARACTERS[sizeof ALLOWED_SPECIAL_CHARACTERS];
}

/**
 * Returns whether the specified character is authorized for an URL path.
 */
bool Url::PathNormalizer::isAllowedCharacter(char c)
{
	return isalnum(c) or isAllowedSpecial(c);
}

/**
 * Stream insertion operator for PathNormalizer: Writes a file path as an URL path.
 */
std::ostream& operator<<(std::ostream& stream, const Url::PathNormalizer& normalizer)
{
	for(size_t i = 0; i < normalizer.path.size(); ++i)
	{
		const char c = normalizer.path[i];
		if(Url::PathNormalizer::isAllowedCharacter(c))
		{
			stream << c;
		}
		else
		{
			stream << '%' << std::hex << static_cast<int>(c) << std::dec;
		}
	}
	return stream;
}

/**
 * Stream extraction operator for PathNormalizer: Reads a file path from an URL path.
 */
std::istream& operator>>(std::istream& stream, const Url::PathNormalizer& normalizer)
{
	char c = stream.get();
	while(Url::PathNormalizer::isAllowedCharacter(c) or c == '%')
	{
		if(c == '%')
		{
			char ordinal[3] = {static_cast<char>(stream.get()), static_cast<char>(stream.get()), '\0'};
			normalizer.path += static_cast<char>(strtol(ordinal, 0, 16));
		}
		else
		{
			normalizer.path += c;
		}
		c = stream.get();
	}
	return stream;
}



///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Url::ParseError class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor. Initializes exception informations.
 * @param url the erroneous URL.
 * @param errorCursor the position of the character where parsing failed.
 */
Url::ParseError::ParseError(const std::string& url, int errorCursor)
	: url(url), errorCursor(errorCursor)
{

}

/**
 * Destructor.
 */
Url::ParseError::~ParseError() throw()
{
}

/**
 * Getter of the "url" property.
 * @return Returns the erroneous URL.
 */
const std::string& Url::ParseError::getUrl() const
{
	return url;
}

/**
 * Getter of the "errorCursor" property.
 * @return Returns the position of the character where parsing failed.
 */
int Url::ParseError::getErrorCursor() const
{
	return errorCursor;
}

/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string Url::ParseError::formatMessage() const
{
	std::ostringstream message;
	message << "The URL: \"" << url << "\" is invalid. There is an error at character "
		<< errorCursor << ": \"" << url.substr(errorCursor) << "\".";
	return message.str();
}


///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::UrlParser class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param url a pointer on the URL object which will take parsing results.
 */
UrlParser::UrlParser(Url* url)
	: url(url)
{
	valid.set();
}

/**
 * Semantic action triggered when a scheme has been successfully parsed.
 * @param begin an iterator to the token first character.
 * @param end an iterator one unit beyond the token last character.
 */
void UrlParser::onScheme(const char* begin, const char* end)
{
	url->scheme.assign(begin, end);
	valid.set(SCHEME, 0);
}

/**
 * Semantic action triggered when a user name has been successfully parsed.
 * @param begin an iterator to the token first character.
 * @param end an iterator one unit beyond the token last character.
 */
void UrlParser::onUser(const char* begin, const char* end)
{
	url->user.assign(begin, end);
	valid.set(AUTH, 0);
}

/**
 * Semantic action triggered when a password has been successfully parsed.
 * @param begin an iterator to the token first character.
 * @param end an iterator one unit beyond the token last character.
 */
void UrlParser::onPassword(const char* begin, const char* end)
{
	url->password.assign(begin, end);
	valid.set(AUTH, 0);
}

/**
 * Semantic action triggered when a host has been successfully parsed.
 * @param begin an iterator to the token first character.
 * @param end an iterator one unit beyond the token last character.
 */
void UrlParser::onHost(const char* begin, const char* end)
{
	url->host.assign(begin, end);
	valid.set(HOSTPORT, 0);
}

/**
 * Semantic action triggered when a port has been successfully parsed.
 * @param value the token value.
 */
void UrlParser::onPort(int value)
{
	url->port = value;
	valid.set(HOSTPORT, 0);
}

/**
 * Semantic action triggered when a path has been successfully parsed.
 * @param begin an iterator to the token first character.
 * @param end an iterator one unit beyond the token last character.
 */
void UrlParser::onPath(const char* begin, const char* end)
{
	url->path.assign(begin, end);
	valid.set(PATH, 0);
}

/**
 * Semantic action triggered when a query item has been successfully parsed.
 * The iterator range encompass the whole item, so it delimits a string of the
 * form key=value.
 * @param begin an iterator to the token first character.
 * @param end an iterator one unit beyond the token last character.
 */
void UrlParser::onQueryItem(const char* begin, const char* end)
{
	std::string queryPair(begin, end);
	int pos = queryPair.find_first_of('=');

	std::string key = queryPair.substr(0, pos);
	std::string value = queryPair.substr(pos + 1);
	url->addQueryItem(key, value);

	valid.set(QUERY, 0);
}

/**
 * Semantic action triggered when possible ambiguity about the previously
 * parsed scheme token is cleared up.
 * Indeed, it may be confused with a relative url path where the first
 * segment would have a colon (:).
 */
void UrlParser::validateScheme(const char*, const char*)
{
	valid.set(SCHEME);
}

/**
 * Semantic action triggered when possible ambiguity about the previously
 * parsed authentication informations is cleared up.
 * Indeed they may be confused with host:port before finding the '@'.
 */
void UrlParser::validateAuth(const char*, const char*)
{
	valid.set(AUTH);
}

/**
 * Semantic action triggered when possible ambiguity about the previously
 * parsed host and port tokens is cleared up.
 */
void UrlParser::validateHostport(const char*, const char*)
{
	valid.set(HOSTPORT);
}

/**
 * Semantic action triggered when possible ambiguity about the previously
 * parsed login informations is cleared up.
 */
void UrlParser::validateLogin(const char*, const char*)
{
	valid.set(LOGIN);
}

/**
 * Semantic action triggered when possible ambiguity about the previously
 * parsed path token is cleared up.
 */
void UrlParser::validatePath(const char*, const char*)
{
	valid.set(PATH);
}

/**
 * Semantic action triggered when possible ambiguity about the previously
 * parsed query is cleared up.
 */
void UrlParser::validateQuery(const char*, const char*)
{
	valid.set(QUERY);
}

/**
 * Semantic action triggered when the whole URL has been successfully parsed.
 * It erases non validated properties.
 */
void UrlParser::validateUrl(const char*, const char*)
{
	if(not valid.test(SCHEME))
	{
		url->scheme.clear();
	}
	if(not valid.test(AUTH))
	{
		url->user.clear();
		url->password.clear();
	}
	if(not valid.test(HOSTPORT))
	{
		url->host.clear();
		url->port = 0;
	}
	if(not valid.test(LOGIN))
	{
		url->user.clear();
		url->password.clear();
		url->host.clear();
		url->port = 0;
	}
	if(not valid.test(PATH))
	{
		url->path.clear();
	}
	if(not valid.test(QUERY))
	{
		url->query.clear();
	}
	valid.set();
}

}

}
