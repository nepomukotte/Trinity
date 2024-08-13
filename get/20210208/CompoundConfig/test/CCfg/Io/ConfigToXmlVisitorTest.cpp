// =============================================================================
// ConfigToXmlVisitorTest.cpp
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
#include "CCfg/Io/ConfigToXmlVisitor.h"
#include "CCfg/Io/IncludeHandler.h"
#include <boost/test/auto_unit_test.hpp>
#include <sstream>
#include <iostream>

namespace CCfg
{
namespace Io
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config2xml_test_data : IoTestConfig
{

	config2xml_test_data()
	{
	}

	~config2xml_test_data()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Io_ConfigToXmlVisitor, config2xml_test_data)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 *
 */
BOOST_AUTO_TEST_CASE(testConfigToXmlVisitor1)
{
	std::ostringstream stream;
	ConfigToXmlVisitor visitor(NoOpIncludeHandler::getInstance(), stream);
	daq.accept(visitor);

	BOOST_CHECK_EQUAL(stream.str(),
		"<Daq>\n"
		"\t<!--The data acquistion duration for a frame-->\n"
		"\t<frameDuration unit=\"ms\">100\n"
		"\t\t<!--10 per second-->\n"
		"\t</frameDuration>\n"
		"\t<readOutFreq unit=\"Mb/s\" range=\"[10, 80]\">25</readOutFreq>\n"
		"\t<Channel id=\"*\">\n"
		"\t\t<isActive>true</isActive>\n"
		"\t\t<!--Data flow parameters for each type of data-->\n"
		"\t\t<Flow id=\"SPE\">\n"
		"\t\t\t<eventRate unit=\"kHz\" range=\"{2, 100}\">100</eventRate>\n"
		"\t\t\t<eventSize unit=\"Bytes\" range=\"{6, 250}\">6</eventSize>\n"
		"\t\t\t<!--Single Photo Electron-->\n"
		"\t\t</Flow>\n"
		"\t\t<Flow id=\"WFA\">\n"
		"\t\t\t<eventRate>2</eventRate>\n"
		"\t\t\t<eventSize>250</eventSize>\n"
		"\t\t</Flow>\n"
		"\t</Channel>\n"
		"\t<Channel id=\"0\" />\n"
		"\t<Channel id=\"1\">\n"
		"\t\t<isActive>false</isActive>\n"
		"\t</Channel>\n"
		"\t<Channel id=\"2\" />\n"
		"</Daq>\n");
}

/**
 *
 */
BOOST_AUTO_TEST_CASE(testConfigToXmlVisitor2)
{
	std::ostringstream stream;
	ConfigToXmlVisitor visitor(NoOpIncludeHandler::getInstance(), stream);
	config.accept(visitor);

	BOOST_CHECK_EQUAL(stream.str(),
		"<Conf>\n"
		"\t<scalar unit=\"Tabs\tand even newlines\nin this unit!!\">"
		"Entities: &lt;&gt;&amp;&apos;&quot;. Other characters:&#x7;&#x1;&#xf;</scalar>\n"
		"</Conf>\n");
}

/**
 *
 */
BOOST_AUTO_TEST_CASE(testConfigToXmlVisitor3)
{
	std::ostringstream stream;
	ConfigToXmlVisitor visitor(NoOpIncludeHandler::getInstance(), stream);
	ambiguityTest.accept(visitor);
	BOOST_CHECK_EQUAL(stream.str(),
		"<AmbiguityTest>\n"
		"\t<intString type=\"string\">123</intString>\n"
		"\t<hexString type=\"string\">0xABC</hexString>\n"
		"\t<floatString type=\"string\">5.0</floatString>\n"
		"\t<boolString type=\"string\">true</boolString>\n"
		"\t<emptyString type=\"string\"></emptyString>\n"
		"\t<emptyObject />\n"
		"\t<intReal>5.0</intReal>\n"
		"</AmbiguityTest>\n");
}

BOOST_AUTO_TEST_SUITE_END()
}
}
