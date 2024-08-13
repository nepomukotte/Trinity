// =============================================================================
// UrlTest.cpp
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

#include "CCfg/Io/Url.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{
namespace Io
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct url_parser_test
{
	url_parser_test()
	{
	}
	~url_parser_test()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Io_Url, url_parser_test)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(testUrl1)
{
	Url url("/home/anonymous/foo.c");
	BOOST_CHECK_MESSAGE(not url.isRelative(), "relative");
	BOOST_CHECK_MESSAGE(not url.hasUser(), "user");
	BOOST_CHECK_MESSAGE(not url.hasPassword(), "passwd");
	BOOST_CHECK_MESSAGE(not url.hasHost(), "host");
	BOOST_CHECK_EQUAL(url.getPort(), 0);
	BOOST_CHECK_EQUAL(url.getPath(), "/home/anonymous/foo.c");
	BOOST_CHECK_MESSAGE(not url.hasQuery(), "query");
	BOOST_CHECK_EQUAL(url.format(), "/home/anonymous/foo.c");

}

BOOST_AUTO_TEST_CASE(testUrl2)
{
	Url url("my-ftp://hello-server/home/blah/foo.c");
	BOOST_CHECK_EQUAL(url.getScheme(), "my-ftp");
	BOOST_CHECK_MESSAGE(not url.hasUser(), "user");
	BOOST_CHECK_MESSAGE(not url.hasPassword(), "passwd");
	BOOST_CHECK_EQUAL(url.getHost(), "hello-server");
	BOOST_CHECK_EQUAL(url.getPort(), 0);
	BOOST_CHECK_EQUAL(url.getPath(), "home/blah/foo.c");
	BOOST_CHECK_MESSAGE(not url.hasQuery(), "query");
	BOOST_CHECK_EQUAL(url.format(), "my-ftp://hello-server/home/blah/foo.c");
}

BOOST_AUTO_TEST_CASE(testUrl3)
{
	Url url("http://user-512:hello720$@hyper-server:12345/user-512/foo.c");
	BOOST_CHECK_EQUAL(url.getScheme(), "http");
	BOOST_CHECK_EQUAL(url.getUser(), "user-512");
	BOOST_CHECK_EQUAL(url.getPassword(), "hello720$");
	BOOST_CHECK_EQUAL(url.getHost(), "hyper-server");
	BOOST_CHECK_EQUAL(url.getPort(), 12345);
	BOOST_CHECK_EQUAL(url.getPath(), "user-512/foo.c");
	BOOST_CHECK_MESSAGE(not url.hasQuery(), "query");
	BOOST_CHECK_EQUAL(url.format(), "http://user-512:hello720$@hyper-server:12345/user-512/foo.c");
}

BOOST_AUTO_TEST_CASE(testUrl4)
{
	Url url("http://ano-nymous91@serv-18.corp-5.com:50000/ano-nymous/test:%201024@85/doc.txt");
	BOOST_CHECK_EQUAL(url.getScheme(), "http");
	BOOST_CHECK_EQUAL(url.getUser(), "ano-nymous91");
	BOOST_CHECK_MESSAGE(not url.hasPassword(), "passwd");
	BOOST_CHECK_EQUAL(url.getHost(), "serv-18.corp-5.com");
	BOOST_CHECK_EQUAL(url.getPort(), 50000);
	BOOST_CHECK_EQUAL(url.getPath(), "ano-nymous/test:%201024@85/doc.txt");
	BOOST_CHECK_MESSAGE(not url.hasQuery(), "passwd");
	BOOST_CHECK_EQUAL(url.format(), "http://ano-nymous91@serv-18.corp-5.com:50000/ano-nymous/test:%201024@85/doc.txt");
}

BOOST_AUTO_TEST_CASE(testUrl5)
{
	Url url("http://192.168.0.1:8080/index.php?search=foo&type=bar");
	BOOST_CHECK_EQUAL(url.getScheme(), "http");
	BOOST_CHECK_MESSAGE(not url.hasUser(), "user");
	BOOST_CHECK_MESSAGE(not url.hasPassword(), "passwd");
	BOOST_CHECK_EQUAL(url.getHost(), "192.168.0.1");
	BOOST_CHECK_EQUAL(url.getPort(), 8080);
	BOOST_CHECK_EQUAL(url.getPath(), "index.php");
	BOOST_CHECK_EQUAL(url.getQueryItem("search"), "foo");
	BOOST_CHECK_EQUAL(url.getQueryItem("type"), "bar");
	BOOST_CHECK_EQUAL(url.format(), "http://192.168.0.1:8080/index.php?search=foo&type=bar");
}

BOOST_AUTO_TEST_CASE(testUrl6)
{
	Url url("test:42-$/result@/zorglub.txt");
	BOOST_CHECK_MESSAGE(url.isRelative(), "relative");
	BOOST_CHECK_MESSAGE(not url.hasUser(), "user");
	BOOST_CHECK_MESSAGE(not url.hasPassword(), "passwd");
	BOOST_CHECK_MESSAGE(not url.hasHost(), "host");
	BOOST_CHECK_EQUAL(url.getPort(), 0);
	BOOST_CHECK_EQUAL(url.getPath(), "test:42-$/result@/zorglub.txt");
	BOOST_CHECK_MESSAGE(not url.hasQuery(), "query");
	BOOST_CHECK_EQUAL(url.format(), "test:42-$/result@/zorglub.txt");
}

BOOST_AUTO_TEST_CASE(testUrl7)
{
	try
	{
		Url url("http://hello-:80/doc.html");
		BOOST_FAIL("Should throw ParseError");
	}
	catch(Url::ParseError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The URL: \"http://hello-:80/doc.html\" is invalid. "
			"There is an error at character 12: \"-:80/doc.html\".");
	}
}

BOOST_AUTO_TEST_CASE(testUrl8)
{
	Url url = Url::fromPath("folder<#>:%/foo bar/$$hello-world!");
	BOOST_CHECK_EQUAL(url.getPath(), "folder%3c%23%3e:%25/foo%20bar/$$hello-world!");
	std::string path = url.toPath();
	BOOST_CHECK_EQUAL(path, "folder<#>:%/foo bar/$$hello-world!");
}

BOOST_AUTO_TEST_SUITE_END()
}
}
