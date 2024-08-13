// =============================================================================
// DocumentTest.cpp
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
#include "CCfg/Io/FileStorage.h"
#include "CCfg/Io/XmlFormat.h"
#include "CCfg/Io/Document.h"
#include "CCfg/Io/Url.h"
#include "CCfg/Object.h"
#include "CCfg/Scalar.h"
#include <boost/test/auto_unit_test.hpp>

#include <sys/stat.h>
#ifdef WIN32
#include <io.h>
#endif
#include <unistd.h>

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
bool compareFiles(const std::string& file1, const std::string& file2);
void copyFile(const std::string& dest, const std::string& src);
void createDir(const char* dir);
int changeDir(const char* dir);

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct document_test : IoTestConfig
{
	document_test()
	{
	}

	~document_test()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Io_Document, document_test)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Test successful load(storage, format) function without includes.
 */
BOOST_AUTO_TEST_CASE(testLoad)
{
	FileStorage file(DATA_DIR "daq.xml", FileStorage::Read);

	Document doc;
	XmlFormat xml;
	Object& config = doc.load(file, xml);
	config.setIncludeUrl(""); // clears dynamically added root url.
	BOOST_CHECK(daq.equals(config));
	BOOST_CHECK_EQUAL(doc.getName(), "daq.xml");
}

/**
 * Test modify and save(destStorage, format)
 */
BOOST_AUTO_TEST_CASE(testSave1)
{
	// open a file
	FileStorage file(DATA_DIR "daq.xml", FileStorage::Read);
	Document doc;
	XmlFormat xml;
	Object& config = doc.load(file, xml);

	// modify it
	config.append(OBJECT_INST, "Channel", "3");
	config("Channel", "3").append(OBJECT_INST, "Flow", "WFA");
	Attribute& eventRate = config("Channel", "3")("Flow", "WFA").appendScalar("eventRate");
	eventRate = 2.4;
	eventRate.setUnit("kHz");

	// save it as a new file
	FileStorage outFile(TEMP_DIR "test_doc4.xml", FileStorage::Write | FileStorage::Truncate);
	doc.save(outFile, xml);
	outFile.close();

	// compare the produced file with a benchmark
	BOOST_CHECK(compareFiles(TEMP_DIR "test_doc4.xml", DATA_DIR "daq_mod.xml"));
}

/**
 * Test modify and save (without include).
 */
BOOST_AUTO_TEST_CASE(testSave2)
{
	// copy the daq1 file
	copyFile(TEMP_DIR "test_doc5.xml", DATA_DIR "daq.xml");

	// open the temp file
	FileStorage file(TEMP_DIR "test_doc5.xml", FileStorage::Read);
	Document doc;
	XmlFormat xml;
	Object& config = doc.load(file, xml);

	// modify it
	config.append(OBJECT_INST, "Channel", "3");
	config("Channel", "3").append(OBJECT_INST, "Flow", "WFA");
	Attribute& eventRate = config("Channel", "3")("Flow", "WFA").appendScalar("eventRate");
	eventRate = 2.4;
	eventRate.setUnit("kHz");

	// save it
	doc.save();

	// compare it with a benchmark
	BOOST_CHECK(compareFiles(TEMP_DIR "test_doc5.xml", DATA_DIR "daq_mod.xml"));
}

/**
 * Test load with include (relative, 1 level nesting) and save flat.
 */
BOOST_AUTO_TEST_CASE(testLoadIncludeSaveFlat)
{
	// copy the daq1 file
	copyFile(TEMP_DIR "test_doc6.xml", DATA_DIR "daq_inc.xml");
	createDir(TEMP_DIR "channels");
	copyFile(TEMP_DIR "channels/default_channel.xml", DATA_DIR "channels/default_channel.xml");

	// load the temp file
	Document doc;
	doc.loadFromFile(TEMP_DIR "test_doc6.xml");
	BOOST_CHECK_EQUAL(doc.getName(), "test_doc6.xml");
	BOOST_CHECK_EQUAL(doc.getUrl(), TEMP_DIR "test_doc6.xml");

	// save it
	{
		FileStorage outFile(TEMP_DIR "test_doc6_full.xml", FileStorage::Write | FileStorage::Truncate);
		XmlFormat xml;
		doc.save(outFile, xml);
		BOOST_CHECK_EQUAL(doc.getName(), "test_doc6_full.xml");
		BOOST_CHECK_EQUAL(doc.getUrl(), TEMP_DIR "test_doc6_full.xml");
	}

	// compare it with a benchmark
	BOOST_CHECK(compareFiles(TEMP_DIR "test_doc6_full.xml", DATA_DIR "daq.xml"));
}

BOOST_AUTO_TEST_CASE(testLoadIncludeCurrentDir)
{
	Document doc2(DATA_DIR "daq_inc.xml");
	BOOST_REQUIRE(changeDir("data") == 0);
	try
	{
		Document doc1("daq_inc.xml");
		BOOST_CHECK(doc1.getConfig()->equals(*doc2.getConfig(), EQUIVALENT));
	}
	catch(...)
	{
		BOOST_REQUIRE(changeDir("..") == 0);
		throw;
	}
	BOOST_REQUIRE(changeDir("..") == 0);
}

/**
 * Test load with include (relative, 1 level nesting) and saving to original files.
 */
BOOST_AUTO_TEST_CASE(testLoadIncludeSave)
{
	// copy the daq file
	copyFile(TEMP_DIR "test_doc7.xml", DATA_DIR "daq_inc.xml");
	createDir(TEMP_DIR "channels");
	copyFile(TEMP_DIR "channels/default_channel.xml", DATA_DIR "channels/default_channel.xml");

	// load the temp file
	Document doc;
	Object& config = doc.loadFromFile(TEMP_DIR "test_doc7.xml");
	BOOST_CHECK_EQUAL(doc.getName(), "test_doc7.xml");
	BOOST_CHECK_EQUAL(doc.getUrl(), TEMP_DIR "test_doc7.xml");

	// modify the base part
	config.append(OBJECT_INST, "Channel", "3");
	config("Channel", "3").append(OBJECT_INST, "Flow", "WFA");
	Attribute& eventRate = config("Channel", "3")("Flow", "WFA").appendScalar("eventRate");
	eventRate = 2.4;
	eventRate.setUnit("kHz");
	// modify the included part
	config("Channel", "*")("Flow", "WFA")("eventSize") = 256;

	// save the config to its original files
	doc.save();
	BOOST_CHECK_EQUAL(doc.getName(), "test_doc7.xml");
	BOOST_CHECK_EQUAL(doc.getUrl(), TEMP_DIR "test_doc7.xml");

	// compare them with a benchmark
	BOOST_CHECK(compareFiles(TEMP_DIR "test_doc7.xml", DATA_DIR "daq_inc_mod.xml"));
	BOOST_CHECK(compareFiles(TEMP_DIR "channels/default_channel.xml", DATA_DIR "default_channel_mod.xml"));
}

/**
 * Test specifying an inexistent storage.
 */
BOOST_AUTO_TEST_CASE(testInexistentStorage)
{
	try
	{
		FileStorage file;
		file.setFileName(DATA_DIR "daq_inexistent.xml");
		file.setOpenMode(FileStorage::Read);

		Document doc;
		XmlFormat xml;
		doc.load(file, xml);
		BOOST_FAIL("Should throw Storage::OpenError");
	}
	catch(Document::IoError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"Error loading the document located at: \"" DATA_DIR "daq_inexistent.xml\".\n"
			"Cannot open the file \"" DATA_DIR "daq_inexistent.xml\" using mode: 1");
	}
}

/**
 * Test loading an ill formed file.
 */
BOOST_AUTO_TEST_CASE(testLoadIllFormed)
{
	try
	{
		FileStorage file;
		file.open(DATA_DIR "daq_ill.xml", FileStorage::Read);

		Document doc;
		XmlFormat xml;
		Object& config = doc.load(file, xml);
		BOOST_CHECK(daq.equals(config));
		BOOST_FAIL("Should throw Format::SyntaxError");
	}
	catch(Document::IoError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"Error loading the document located at: \"" DATA_DIR "daq_ill.xml\".\n"
			"The text stream is not valid XML.\n"
			"At line 3 column 5: Error reading end tag. (error #9).");
	}
}

/**
 * Test unhandled URL.
 */
BOOST_AUTO_TEST_CASE(testUnhandledUrl)
{
	try
	{
		Document doc;
		doc.loadFromUrl(Url("file://host/doc.xml"));
		BOOST_FAIL("Should throw IoError");
	}
	catch(Document::IoError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"Error loading the document located at: \"file://host/doc.xml\".\n"
			"Cannot build a Format from the URL: \"file://host/doc.xml\".");
	}
}

/**
 * Test recursive include
 */
BOOST_AUTO_TEST_CASE(testRecursiveInclude)
{
	try
	{
		Document doc;
		doc.loadFromFile(DATA_DIR "daq_rec.xml");
		BOOST_FAIL("Should throw IoError");
	}
	catch(Document::IoError& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"Error loading the document located at: \"" DATA_DIR "daq_rec.xml\".\n"
			"The file \"" DATA_DIR  "default_channel_rec.xml\" cannot include the file "
			"\"" DATA_DIR "daq_rec.xml\" because this would begin an infinite recursion.");
	}
}

/**
 * Test include context url computation
 */
BOOST_AUTO_TEST_CASE(testIncludeUrl)
{
	Object config("Daq");
	Object* a = new Object("A");
	Object* b = new Object("B");
	Object* c = new Object("C");
	Object* d = new Object("D");
	Object* e = new Object("E");
	Scalar* x = new Scalar("x");
	config.append(a);
	a->append(b);
	b->append(c);
	c->append(d);
	d->append(e);
	e->append(x);
	config.setIncludeUrl("/home/me/foo.xml");
	a->setIncludeUrl("crates/default.xml");
	c->setIncludeUrl("boards/myboard.xml");
	d->setIncludeUrl("/home/me/global/memory.xml");
	e->setIncludeUrl("banks/default.xml");

	Url url = Document::computeQualifiedIncludeUrl(config);
	BOOST_CHECK_MESSAGE(url.isSameResource(Url::fromPath("/home/me/foo.xml")), "1");

	url = Document::computeQualifiedIncludeUrl(*a);
	BOOST_CHECK_MESSAGE(url.isSameResource(Url::fromPath("/home/me/crates/default.xml")), "2");

	url = Document::computeQualifiedIncludeUrl(*b);
	BOOST_CHECK_MESSAGE(url.isSameResource(Url()), "3");

	url = Document::computeQualifiedIncludeUrl(*c);
	BOOST_CHECK_MESSAGE(url.isSameResource(Url::fromPath("/home/me/crates/boards/myboard.xml")), "4");

	url = Document::computeQualifiedIncludeUrl(*d);
	BOOST_CHECK_MESSAGE(url.isSameResource(Url::fromPath("/home/me/global/memory.xml")), "5");

	url = Document::computeQualifiedIncludeUrl(*e);
	BOOST_CHECK_MESSAGE(url.isSameResource(Url::fromPath("/home/me/global/banks/default.xml")), "6");
}

BOOST_AUTO_TEST_CASE(testDetach)
{
	Document doc;
	Object* configPtr = new Object("Daq");
	doc.setConfig(std::auto_ptr<Object>(configPtr));
	std::auto_ptr<Object> config = doc.detachConfig();
	BOOST_CHECK_EQUAL(doc.getConfig(), static_cast<Object*>(0));
	BOOST_CHECK_EQUAL(config.get(), configPtr);
}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////////////////////////////////////////////
// Utility functions
///////////////////////////////////////////////////////////////////////////////

bool compareFiles(const std::string& file1, const std::string& file2)
{
	FileStorage tempFile(file1, FileStorage::Read);
	FileStorage testFile(file2, FileStorage::Read);

	std::string str1, str2;
	tempFile.readAll(str1);
	testFile.readAll(str2);
	return str1 == str2;
}

void copyFile(const std::string& dest, const std::string& src)
{
	FileStorage inFile(src, FileStorage::Read);
	std::string str;
	inFile.readAll(str);
	FileStorage outFile(dest, FileStorage::Write | FileStorage::Truncate);
	outFile.write(str);
}

void createDir(const char* dir)
{
#if defined(VXWORKS) || defined(WIN32)
	mkdir(dir);
#else
	mkdir(dir, S_IRWXU);
#endif
}

int changeDir(const char* dir)
{
	return chdir(dir);
}

}
}
