// =============================================================================
// XmlFormatTest.cpp
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

#include "Config.h"
#include "CCfg/Io/XmlFormat.h"
#include "CCfg/Io/FileStorage.h"
#include "CCfg/Io/BufferStorage.h"
#include "CCfg/Io/IncludeHandler.h"
#include "CCfg/Object.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{
namespace Io
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct xml_format_test : IoTestConfig
{

	xml_format_test()
	{
	}

	~xml_format_test()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Io_XmlFormat, xml_format_test)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check whether parsing the formatted config in a file returns the original.
 */
BOOST_AUTO_TEST_CASE(testXmlFormat1)
{
	BufferStorage buffer;
	XmlFormat xml;
	xml.format(&daq, buffer, NoOpIncludeHandler::getInstance());

	std::auto_ptr<Object> daq2 = xml.parse(buffer, NoOpIncludeHandler::getInstance());

	BOOST_CHECK(daq.equals(*daq2));
}

/**
 * Check whether parsing an invalid XML stream issues a SyntaxError exception.
 */
BOOST_AUTO_TEST_CASE(testXmlFormat2)
{
	try
	{
		XmlFormat xml;
		BufferStorage buffer;
		buffer.write(
			"<Daq>\n"
			"\t<frameDuration unit=\"ms\">100</frame>\n"
			"</Daq>\n");
		std::auto_ptr<Object> daq2 = xml.parse(buffer, NoOpIncludeHandler::getInstance());
		BOOST_FAIL("Should throw SyntaxError");
	}
	catch(XmlFormat::SyntaxError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The text stream is not valid XML.\n"
			"At line 2 column 33: Error reading end tag. (error #9).");
	}
}

BOOST_AUTO_TEST_SUITE_END()
}
}
