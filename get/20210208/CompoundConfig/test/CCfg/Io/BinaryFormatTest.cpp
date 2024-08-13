// =============================================================================
// BinaryFormatTest.cpp
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Patrick Sizun (patrick.sizun@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "Config.h"
#include "CCfg/Io/BinaryFormat.h"
#include "CCfg/Io/FileStorage.h"
#include "CCfg/Io/BufferStorage.h"
#include "CCfg/Io/IncludeHandler.h"
#include "CCfg/Object.h"
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

struct binary_format_test : IoTestConfig
{

	binary_format_test()
	{
	}

	~binary_format_test()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Io_BinaryFormat, binary_format_test)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check whether formatting the configuration into a file
 *  and parsing this file returns the original (except for comments).
 */
BOOST_AUTO_TEST_CASE(testBinaryFormat1)
{
	FileStorage out(TEMP_DIR "daq.bcfg",
			FileStorage::Write | FileStorage::Truncate | FileStorage::Binary);
	BinaryFormat binary;
	binary.format(&daq, out, NoOpIncludeHandler::getInstance());

	BinaryFormat binary2;
	FileStorage in(TEMP_DIR "daq.bcfg", FileStorage::Read | FileStorage::Binary);
	std::auto_ptr<Object> daq2(binary2.parse(in, NoOpIncludeHandler::getInstance()));

	BOOST_CHECK(daq.equals(*daq2, EQUIVALENT));
}

/**
 * Check whether formatting the configuration into a buffer
 *  and parsing this buffer's contents returns the original (except for comments).
 */
BOOST_AUTO_TEST_CASE(testBinaryFormat2)
{
	BufferStorage buffer;
	BinaryFormat binary;
	binary.format(&daq, buffer, NoOpIncludeHandler::getInstance());

	BinaryFormat binary2;
	std::auto_ptr<Object> daq2 = binary2.parse(buffer, NoOpIncludeHandler::getInstance());

	BOOST_CHECK(daq.equals(*daq2, EQUIVALENT));
}

/**
 * Check whether configuration is identical to configuration from
 * pre-existing binary configuration file (except for comments).
 */
BOOST_AUTO_TEST_CASE(testBinaryFormat3)
{
	BinaryFormat binary;
	FileStorage in(DATA_DIR "daq.bcfg", FileStorage::Read | FileStorage::Binary);
	std::auto_ptr<Object> daq2(binary.parse(in, NoOpIncludeHandler::getInstance()));

	BOOST_CHECK(daq.equals(*daq2, EQUIVALENT));
}

BOOST_AUTO_TEST_SUITE_END()
}
}
