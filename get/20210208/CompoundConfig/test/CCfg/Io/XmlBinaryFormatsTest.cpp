// =============================================================================
// XmlBinaryFormatsTest.cpp
// -----------------------------------------------------------------------------
//  $Id: XmlBinaryFormatsTest.cpp 1170 2013-11-21 13:26:34Z fchateau $
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Patrick Sizun (patrick.sizun@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "Config.h"
#include "CCfg/Io/FileStorage.h"
#include "CCfg/Io/BinaryFormat.h"
#include "CCfg/Io/XmlFormat.h"
#include "CCfg/Io/IncludeHandler.h"
#include "CCfg/Io/Url.h"
#include "CCfg/Object.h"
#include "CCfg/Scalar.h"
#include <boost/test/auto_unit_test.hpp>

#ifndef VXWORKS
#define TEMP_DIR "/tmp/"
#else
#define TEMP_DIR "tmp/"
#endif

#define DATA_DIR "data/"

namespace CCfg
{
namespace Io
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct xml_binary_formats_test : IoTestConfig
{
	xml_binary_formats_test()
	{
	}

	~xml_binary_formats_test()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Io_XmlBinaryFormat, xml_binary_formats_test)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Test successful load(storage, format) function without includes.
 */
BOOST_AUTO_TEST_CASE(testXmlBinaryFormat1)
{
	FileStorage binFile(DATA_DIR "daq2.bcfg", FileStorage::Read | FileStorage::Binary);
	BinaryFormat binary;
	std::auto_ptr<Object> binConfig(binary.parse(binFile, NoOpIncludeHandler::getInstance()));

	FileStorage xmlFile(DATA_DIR "daq2.ccfg", FileStorage::Read );
	XmlFormat xml;
	std::auto_ptr<Object> xmlConfig(xml.parse(xmlFile, NoOpIncludeHandler::getInstance()));
/*
	FileStorage outFile(TEMP_DIR "daq2.bcfg", FileStorage::Write | FileStorage::Binary);
	BinaryFormat binary;
	binary.format(xmlConfig.get(), outFile, NoOpIncludeHandler::getInstance());
*/
	BOOST_CHECK(binConfig->equals(*xmlConfig, EQUIVALENT));
}

BOOST_AUTO_TEST_SUITE_END()
}
}
