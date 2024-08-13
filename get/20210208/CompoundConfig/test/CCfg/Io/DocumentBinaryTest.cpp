// =============================================================================
// DocumentBinaryTest.cpp
// -----------------------------------------------------------------------------
//  $Id$
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
#include "CCfg/Io/Document.h"
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

struct document_binary_test : IoTestConfig
{
	document_binary_test()
	{
	}

	~document_binary_test()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Document_w_BinaryFormat, document_binary_test)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Test successful load(storage, format) function without includes.
 */
BOOST_AUTO_TEST_CASE(testLoadBinary)
{
	FileStorage file(DATA_DIR "daq.bcfg", FileStorage::Read | FileStorage::Binary);

	Document doc;
	BinaryFormat binary;
	Object& config = doc.load(file, binary);
	config.setIncludeUrl(""); // clears dynamically added root url.
	BOOST_CHECK(daq.equals(config, EQUIVALENT));
	BOOST_CHECK_EQUAL(doc.getName(), "daq.bcfg");
}

/**
 * Test specifying an inexistent storage.
 */
BOOST_AUTO_TEST_CASE(testInexistentBinaryStorage)
{
	try
	{
		FileStorage file;
		file.setFileName(DATA_DIR "daq_inexistent.bcfg");
		file.setOpenMode(FileStorage::Read | FileStorage::Binary);

		Document doc;
		BinaryFormat binary;
		doc.load(file, binary);
		BOOST_FAIL("Should throw Storage::OpenError");
	}
	catch(Document::IoError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"Error loading the document located at: \"" DATA_DIR "daq_inexistent.bcfg\".\n"
			"Cannot open the file \"" DATA_DIR "daq_inexistent.bcfg\" using mode: 17");
	}
}

/**
 * Test loading an ill formed file.
 */
BOOST_AUTO_TEST_CASE(testIllFormedFile)
{
	try
	{
		FileStorage file;
		file.open(DATA_DIR "daq2_ill.bcfg", FileStorage::Read);

		Document doc;
		BinaryFormat binary;
		doc.load(file, binary);
		BOOST_FAIL("Should throw Format::SemanticError");
	}
	catch(Document::IoError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"Error loading the document located at: \"" DATA_DIR "daq2_ill.bcfg\".\n"
			"The binary document is not a valid configuration.");
	}
}

/**
 * Test loading an empty file.
 */
BOOST_AUTO_TEST_CASE(testLoadEmpty)
{
	try
	{
		FileStorage emptyFile(TEMP_DIR "daq_empty.bcfg",
				FileStorage::Write | FileStorage::Truncate | FileStorage::Binary);
		emptyFile.close();

		FileStorage file;
		file.open(TEMP_DIR "daq_empty.bcfg", FileStorage::Read | FileStorage::Binary);

		Document doc;
		BinaryFormat binary;
		doc.load(file, binary);
		BOOST_FAIL("Should throw Format::SemanticError");
	}
	catch(Document::IoError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"Error loading the document located at: \"" TEMP_DIR "daq_empty.bcfg\".\n"
			"The binary document is not a valid configuration.");
	}
}

BOOST_AUTO_TEST_SUITE_END()
}
}
