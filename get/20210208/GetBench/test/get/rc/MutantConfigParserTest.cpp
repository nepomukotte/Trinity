/*
 * @file MutantConfigParserTest.cpp
 * $Id: MutantConfigParserTest.cpp 1794 2021-02-03 11:41:04Z psizun $
 */

#include "get/rc/MutantConfigParser.h"
#include "CCfg/Object.h"
#include "CCfg/Io/Document.h"
#include "boost/test/auto_unit_test.hpp"

#define DATA_DIR "test/data/"

namespace get
{
namespace rc
{

BOOST_AUTO_TEST_CASE(L1BTimeoutTest)
{
	CCfg::Io::Document doc;
	CCfg::Object & obj = doc.loadFromFile(DATA_DIR "mutant-config.xcfg");
	CCfg::CConfig cfg(&obj);
	cfg("TwoProtonDecayMode")("isEnabled").setBoolValue(true);


	// Check default timeout is 2.56 us
	MutantConfigParser parser;
	parser.parseConfig(cfg);
	BOOST_CHECK_EQUAL(parser.level1bTimeout_2560ns, 1u);

	// Check 0 leads to null timeout
	cfg("Level", 1)("L1B")("timeout").setRealValue(0.0);
	parser.parseConfig(cfg);
	BOOST_CHECK_EQUAL(parser.level1bTimeout_2560ns, 0u);
}

BOOST_AUTO_TEST_CASE(ExternalClockSelection)
{
	CCfg::Io::Document doc;
	CCfg::Object & obj = doc.loadFromFile(DATA_DIR "mutant-config.xcfg");
	CCfg::CConfig cfg(&obj);
	cfg("GlobalMasterClock")("externalClockSelection").setStringValue("VHDCI");


	// Check default external clock is VHDCI
	MutantConfigParser parser;
	parser.parseConfig(cfg);
	BOOST_CHECK_EQUAL(parser.externalGmcClockSelection, 2u);

	// Check FireWire 1 input external clock
	cfg("GlobalMasterClock")("externalClockSelection").setStringValue("FireWireA");
	parser.parseConfig(cfg);
	BOOST_CHECK_EQUAL(parser.externalGmcClockSelection, 0u);

	// Check FireWire 2 input external clock
	cfg("GlobalMasterClock")("externalClockSelection").setStringValue("FireWire2");
	parser.parseConfig(cfg);
	BOOST_CHECK_EQUAL(parser.externalGmcClockSelection, 1u);

	// Check bad input external clock
	cfg("GlobalMasterClock")("externalClockSelection").setStringValue("FireWire3");
	parser.parseConfig(cfg);
	BOOST_CHECK_EQUAL(parser.externalGmcClockSelection, 2u);

	// Check Optical input external clock
	cfg("GlobalMasterClock")("externalClockSelection").setStringValue("Optical");
	parser.parseConfig(cfg);
	BOOST_CHECK_EQUAL(parser.externalGmcClockSelection, 3u);
}


BOOST_AUTO_TEST_CASE(AmcTagOnly)
{
	CCfg::Io::Document doc;
	CCfg::Object & obj = doc.loadFromFile(DATA_DIR "mutant-config.xcfg");
	CCfg::CConfig cfg(&obj);
	cfg("AmcTagOnly")("ATO", 3).setBoolValue(true);

	// Check AMC_TAG_ONLY
	MutantConfigParser parser;
	parser.parseConfig(cfg);
	BOOST_CHECK_EQUAL(parser.amcTagOnly, 4u);

	cfg("AmcTagOnly")("ATO", 10).setBoolValue(true);
	parser.parseConfig(cfg);
	BOOST_CHECK_EQUAL(parser.amcTagOnly, 0x204u);
}

} /* namespace rc */
} /* namespace get */




