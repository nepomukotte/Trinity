// =============================================================================
// XmlToConfigVisitorTest.cpp
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
#include "CCfg/Io/XmlToConfigVisitor.h"
#include "CCfg/Io/ConfigToXmlVisitor.h"
#include "CCfg/Io/IncludeHandler.h"
#include "tinyxml.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{
namespace Io
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct xml2config_test_data : IoTestConfig
{

	xml2config_test_data()
	{
	}

	~xml2config_test_data()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Io_XmlToConfigVisitor, xml2config_test_data)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Checks whether parsing a writen configuration do not lose informations.
 */
BOOST_AUTO_TEST_CASE(testXmlToConfigVisitor1)
{
	// write a configuration
	std::ostringstream stream;
	ConfigToXmlVisitor formattingVisitor(NoOpIncludeHandler::getInstance(), stream);
	daq.accept(formattingVisitor);

	// parse it
	TiXmlDocument doc;
	doc.Parse(stream.str().c_str());
	XmlToConfigVisitor parsingVisitor;
	doc.Accept(&parsingVisitor);

	// ensure we get the same configuration as the original.
	BOOST_CHECK(daq.equals(*parsingVisitor.getBuiltConfig()));
}

/**
 * Checks whether an error is issued when an XML element has both txt and
 * children XML elements.
 */
BOOST_AUTO_TEST_CASE(testXmlToConfigVisitor2)
{
	try
	{
		TiXmlDocument doc;
		doc.Parse(
			"\t<Daq>\n"
			"\t\t<Channel id=\"*\">\n"
			"\t\t\t42\n"
			"\t\t\t<isActive>true</isActive>\n"
			"\t\t</Channel>\n"
			"\t</Daq>\n");
		XmlToConfigVisitor visitor;
		doc.Accept(&visitor);
		BOOST_FAIL("Should throw SemanticError");
	}
	catch(XmlFormat::SemanticError& e)
	{
		BOOST_CHECK_EQUAL(e.getErrorId(), XmlSemError::MIXING_TEXT_AND_ELEMENTS);
		BOOST_CHECK_EQUAL(e.getElement(), "<Channel id=\"*\">");
		BOOST_CHECK_EQUAL(e.getErrorRow(), 2);
		BOOST_CHECK_EQUAL(e.getErrorCol(), 9);
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The Xml document is not a valid configuration.\n"
			"The element: <Channel id=\"*\"> at line 2 column 9 is causing a problem.\n"
			"It cannot have both text and children XML elements.");
	}
}

/**
 * Checks whether an error is issued when an XML element has split text.
 */
BOOST_AUTO_TEST_CASE(testXmlToConfigVisitor3)
{
	try
	{
		TiXmlDocument doc;
		doc.Parse(
			"<Daq>\n"
			"\t<frameDuration unit=\"ms\">1<!--10 per second-->00</frameDuration>\n"
			"</Daq>\n");
		XmlToConfigVisitor visitor(NoOpIncludeHandler::getInstance());
		doc.Accept(&visitor);
		BOOST_FAIL("Should throw SemanticError");
	}
	catch(XmlFormat::SemanticError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The Xml document is not a valid configuration.\n"
			"The element: <frameDuration unit=\"ms\"> "
			"at line 2 column 5 is causing a problem.\n"
			"The text of a scalar must not be split.");
	}
}

/**
 * Check whether an error is issued when an XML element standing for a config
 * scalar has an invalid attribute.
 */
BOOST_AUTO_TEST_CASE(testXmlToConfigVisitor4)
{
	try
	{
		TiXmlDocument doc;
		doc.Parse(
			"<Daq>\n"
			"\t<frameDuration unit=\"ms\" foo=\"bar\">100</frameDuration>\n"
			"</Daq>\n");
		XmlToConfigVisitor visitor(NoOpIncludeHandler::getInstance());
		doc.Accept(&visitor);
		BOOST_FAIL("Should throw SemanticError");
	}
	catch(XmlFormat::SemanticError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The Xml document is not a valid configuration.\n"
			"The element: <frameDuration unit=\"ms\" foo=\"bar\"> "
			"at line 2 column 5 is causing a problem.\n"
			"It has an XML attribute that is not allowed by configuration "
			"specifications.");
	}
}

/**
 * Check whether an error is issued when an XML element standing for a config
 * object has an invalid attribute.
 */
BOOST_AUTO_TEST_CASE(testXmlToConfigVisitor5)
{
	try
	{
		TiXmlDocument doc;
		doc.Parse(
			"<Daq>\n"
			"\t<Channel unit=\"ms\" />\n"
			"</Daq>\n");
		XmlToConfigVisitor visitor(NoOpIncludeHandler::getInstance());
		doc.Accept(&visitor);
		BOOST_FAIL("Should throw SemanticError");
	}
	catch(XmlFormat::SemanticError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The Xml document is not a valid configuration.\n"
			"The element: <Channel unit=\"ms\" /> "
			"at line 2 column 5 is causing a problem.\n"
			"It has an XML attribute that is not allowed by configuration "
			"specifications.");
	}
}

/**
 * Check if special characters are restored correctly
 */
BOOST_AUTO_TEST_CASE(testXmlToConfigVisitor6)
{
	// write a configuration
	std::ostringstream stream;
	ConfigToXmlVisitor formattingVisitor(NoOpIncludeHandler::getInstance(), stream);
	config.accept(formattingVisitor);

	// parse it
	TiXmlDocument doc;
	doc.Parse(stream.str().c_str());
	XmlToConfigVisitor parsingVisitor(NoOpIncludeHandler::getInstance());
	doc.Accept(&parsingVisitor);
	std::auto_ptr<Object> config2 = parsingVisitor.getBuiltConfig();

	// ensure we get the same configuration as the original.
	BOOST_CHECK(config.equals(*config2));
}

/**
 * Check if type specifications are interpreted correctly.
 */
BOOST_AUTO_TEST_CASE(testXmlToConfigVisitor7)
{
	// write a configuration
	std::ostringstream stream;
	ConfigToXmlVisitor visitor(NoOpIncludeHandler::getInstance(), stream);
	ambiguityTest.accept(visitor);

	// parse it
	TiXmlDocument doc;
	doc.Parse(stream.str().c_str());
	XmlToConfigVisitor parsingVisitor(NoOpIncludeHandler::getInstance());
	doc.Accept(&parsingVisitor);
	std::auto_ptr<Object> config2 = parsingVisitor.getBuiltConfig();

	// ensure we get the same configuration as the original.
	BOOST_CHECK(ambiguityTest.equals(*config2));
}

/**
 * Check whether an exception is thrown when the specified scalar type is invalid.
 */
BOOST_AUTO_TEST_CASE(testXmlToConfigVisitor8)
{
	try
	{
		TiXmlDocument doc;
		doc.Parse(
			"<Daq>\n"
			"\t<myvalue type=\"complex\">5+3i</myvalue>\n"
			"</Daq>\n");
		XmlToConfigVisitor visitor(NoOpIncludeHandler::getInstance());
		doc.Accept(&visitor);
		BOOST_FAIL("Should throw SemanticError");
	}
	catch(XmlFormat::SemanticError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The Xml document is not a valid configuration.\n"
			"The element: <myvalue type=\"complex\"> "
			"at line 2 column 5 is causing a problem.\n"
			"The scalar type is invalid.");
	}
}

BOOST_AUTO_TEST_SUITE_END()
}
}
