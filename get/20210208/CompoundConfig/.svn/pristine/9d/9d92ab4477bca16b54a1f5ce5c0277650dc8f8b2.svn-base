// =============================================================================
// AttributeTest.cpp
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

#include "View/Configs.h"
#include "CCfg/CConfig.h"
#include "CCfg/Object.h"
#include "CCfg/Scalar.h"
#include "CCfg/Attribute.h"
#include "CCfg/Range.h"
#include "CCfg/UnsupportedOperation.h"
#include <boost/test/auto_unit_test.hpp>
#include <string.h>

#define CATCH_ATTRIBUTE_NOT_FOUND(path, expr) \
	try	\
	{	\
		expr;	\
		BOOST_FAIL("No exception thrown."); \
	} catch(CConfig::NotFound& e)	\
	{	\
		BOOST_CHECK_EQUAL(e.getConfigPath().toString(), path);	\
		BOOST_CHECK_EQUAL(e.getMessage(), "The config element " path " cannot be resolved to a concrete attribute."); \
	} catch(std::logic_error& e) \
	{	\
		throw;	\
	} catch(...)	\
	{	\
		BOOST_FAIL("Bad exception type.");	\
	}


namespace CCfg
{
using View::Result;

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_cconfig
{
	Config1 data1;
	Config2 data2;
	Config3 data3;
	Config4 data4;
	Config5 data5;

	config_cconfig()
	{

	}

	~config_cconfig()
	{

	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_CConfig, config_cconfig)

///////////////////////////////////////////////////////////////////////////////
// fixture 1 test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check if search with one backtracking works.
 */
BOOST_AUTO_TEST_CASE(testSearch1Back)
{
	CConfig config(&data1.cfg);
	Result res = config("A", "1")("x").resolve();
	BOOST_CHECK_EQUAL(res.getAttribute(), data1.a_x);
	BOOST_CHECK_EQUAL(res.getBacktrackLevel(), 1);
}

/**
 * Check if direct attribute access (no backtracking at all) works.
 */
BOOST_AUTO_TEST_CASE(testSearchNoBack)
{
	CConfig config(&data1.cfg);
	Result res = config("A", "*")("x").resolve();
	BOOST_CHECK_EQUAL(res.getAttribute(), data1.a_x);
	BOOST_CHECK_EQUAL(res.getBacktrackLevel(), 0);
}

/**
 * Check if backtracking beyond the base attribute works.
 */
BOOST_AUTO_TEST_CASE(testSearchFullBack)
{
	CConfig config(data1.a1);
	Attribute* attr = config("x").search();
	BOOST_CHECK_EQUAL(attr, data1.a_x);
}

/**
 * Check if backtracking beyond the base attribute fails when necessary.
 */
BOOST_AUTO_TEST_CASE(testSearchFullBack_Failure)
{
	CConfig config(data1.a1);
	CATCH_ATTRIBUTE_NOT_FOUND("cfg.A[1].w",
		Attribute& attr = config("w").find();
		attr = 1
	);
}

/**
 * Check if searching a non existent attribute works.
 */
BOOST_AUTO_TEST_CASE(testSearchNotFound)
{
	CConfig config(&data1.cfg);
	CATCH_ATTRIBUTE_NOT_FOUND("cfg.A[2].x",
		Attribute& attr = config("A", "2")("x").find();
		attr = 1
	);
}

/**
 * Check if attribute value overriding works.
 */
BOOST_AUTO_TEST_CASE(testSearchOverriding)
{
	CConfig config(&data1.cfg);
	Attribute& attr = config("A", "1")("y").find();
	BOOST_CHECK_EQUAL(attr.getIntValue(), 5);
}

/**
 * Check if executing 2 consecutive searches works.
 */
BOOST_AUTO_TEST_CASE(testSearchConsecutive)
{
	CConfig config(&data1.cfg);
	Attribute* attr = config("A", "1")("x").search();
	BOOST_CHECK_EQUAL(attr, data1.a_x);
	attr = config("A", "1").search();
	BOOST_CHECK_EQUAL(attr, data1.a1);
}

/**
 * Check if the trivial search works.
 */
BOOST_AUTO_TEST_CASE(testSearchDummy)
{
	CConfig config(data1.a1);
	Result res = config.resolve();
	BOOST_CHECK_EQUAL(res.getAttribute(), data1.a1);
	BOOST_CHECK_EQUAL(res.getBacktrackLevel(), 0);
}

/**
 * Check if iteration works
 */
BOOST_AUTO_TEST_CASE(testIterator)
{
	CConfig config(&data1.cfg);
	CConfig::Iterator it = config("A", "1").iterate();

	BOOST_CHECK_MESSAGE(it.hasNext(), "next1");
	Attribute* attr = it.next().search();
	BOOST_CHECK_MESSAGE(attr == data1.a_x, "A[*].x");
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 1);

	BOOST_CHECK_MESSAGE(it.hasNext(), "next2");
	it.next();
	attr = it.get().search();
	BOOST_CHECK_MESSAGE(attr == data1.a1y, "A[1].y");
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 0);

	BOOST_CHECK_MESSAGE(not it.hasNext(), "no more attributes");
}


///////////////////////////////////////////////////////////////////////////////
// fixture 2 test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check if searching an object without backtracking works.
 */
BOOST_AUTO_TEST_CASE(testSearchObjectNoBack)
{
	CConfig config(&data2.a);
	Result res = config("B", "0")("C", "0")("D", "0").resolve();
	BOOST_CHECK_EQUAL(res.getAttribute(), data2.ab0c0d0);
	BOOST_CHECK_EQUAL(res.getBacktrackLevel(), 0);
}

/**
 * Check if searching a scalar without backtracking works.
 */
BOOST_AUTO_TEST_CASE(testSearchScalarNoBack)
{
	CConfig config(&data2.a);
	Attribute* attr = config("B", "0")("C", "0")("D", "0")("z").search();
	BOOST_CHECK_EQUAL(attr, data2.z);
}

/**
 * Check if searching a scalar with a 1 backtracking works.
 */
BOOST_AUTO_TEST_CASE(testSearchScalar1Back)
{
	CConfig config(&data2.a);
	Result res = config("B", "0")("C", "0")("D", "0")("y").resolve();
	BOOST_CHECK_EQUAL(res.getAttribute(), data2.y);
	BOOST_CHECK_EQUAL(res.getBacktrackLevel(), 1);
}

/**
 * Check if searching a scalar with a 3 backtracking works.
 */
BOOST_AUTO_TEST_CASE(testSearchScalar3Back)
{
	CConfig config(&data2.a);
	Result res = config("B", "0")("C", "0")("D", "0")("x").resolve();
	BOOST_CHECK_EQUAL(res.getAttribute(), data2.x);
	BOOST_CHECK_EQUAL(res.getBacktrackLevel(), 3);
}

/**
 * Check if backtracking beyond the base attribute works.
 */
BOOST_AUTO_TEST_CASE(testSearchScalarFullBack)
{
	CConfig config(data2.ab0c0d0);
	Result res = config("x").resolve();
	BOOST_CHECK_EQUAL(res.getAttribute(), data2.x);
	BOOST_CHECK_EQUAL(res.getBacktrackLevel(), 3);
}

/**
 * Check if searching a non existent attribute works.
 */
BOOST_AUTO_TEST_CASE(testSearchNotFound2)
{
	CConfig config(&data2.a);
	CATCH_ATTRIBUTE_NOT_FOUND("A.B[0].C[0].D[0].w",
		Attribute& attr = config("B", "0")("C", "0")("D", "0")("w").find();
		attr = 1
	);
}

/**
 * Check if iterating on a non object fails as expected.
 */
BOOST_AUTO_TEST_CASE(testIterateFailure)
{
	CConfig config(&data2.a);

	try
	{
		CConfig::Iterator it = config("B", "0")("C", "0")("D", "0")("x").iterate();
		it.hasNext();
		BOOST_FAIL("No exception thrown.");
	}
	catch(UnsupportedOperation& e)
	{

	}
	catch(std::logic_error& e)
	{
		throw;
	}
	catch(...)
	{
		BOOST_FAIL("Bad exception type.");
	}

}

/**
 * Check if reseting the base attribute works.
 */
BOOST_AUTO_TEST_CASE(testSetBaseAttribute)
{
	CConfig config(&data2.a);
	config.setBaseAttribute(data2.ab0c0d0);
	Attribute* attr = config("x").search();
	BOOST_CHECK_EQUAL(attr, data2.x);
}

/**
 * Check if stream insertion works.
 */
BOOST_AUTO_TEST_CASE(testOstream)
{
	std::ostringstream str;
	CConfig config(&data2.a);
	str << config("B", "0")("C", "0")("D", "0")("x");
	BOOST_CHECK_EQUAL(str.str(), "4.5");
}

/**
 * Check if getName works
 */
BOOST_AUTO_TEST_CASE(testGetName)
{
	CConfig config(&data2.a);
	BOOST_CHECK_EQUAL(config("B", "0")("C", "0")("D", "0").getName(), "D");
	BOOST_CHECK_EQUAL(config("B", "0")("C", "0")("D", "0")("x").getName(), "x");
}

/**
 * Check if getIndex works
 */
BOOST_AUTO_TEST_CASE(testGetIndex)
{
	CConfig config(&data2.a);
	BOOST_CHECK_EQUAL(config("B", "0")("C", "0")("D", "0").getIndex(), "0");
	BOOST_CHECK_EQUAL(config("B", "0")("C", "0")("D", "0")("x").getIndex(), "");
}

/**
 * Check if getPreComment works
 */
BOOST_AUTO_TEST_CASE(testGetPreComment)
{
	CConfig config(data2.ab0c0d0);
	BOOST_CHECK_EQUAL(config("x").getPreComment().getText(), "hello");
}

/**
 * Check if getPostComment works
 */
BOOST_AUTO_TEST_CASE(testGetPostComment)
{
	CConfig config(data2.ab0c0d0);
	BOOST_CHECK_EQUAL(config("x").getPostComment().getText(), "world");
}

/**
 * Check if getIncludeUrl works
 */
BOOST_AUTO_TEST_CASE(testGetIncludeUrl)
{
	CConfig config(data2.ab0c0);
	BOOST_CHECK_EQUAL(config("D", "0").getIncludeUrl(), "foo.xml");
}

/**
 * Check if setPreComment works
 */
BOOST_AUTO_TEST_CASE(testGetText)
{
	CConfig config(data2.ab0c0d0);
	config("x").setPreComment(Comment("foo"));
	BOOST_CHECK_EQUAL(data2.x->getPreComment().getText(), "foo");
}

/**
 * Check if setPostComment works
 */
BOOST_AUTO_TEST_CASE(testSetPostComment)
{
	CConfig config(data2.ab0c0d0);
	config("x").setPostComment(Comment("bar"));
	BOOST_CHECK_EQUAL(data2.x->getPostComment().getText(), "bar");
}

/**
 * Check if setIncludeUrl works
 */
BOOST_AUTO_TEST_CASE(testSetIncludeUrl)
{
	CConfig config(data2.ab0c0);
	config("D", 0).setIncludeUrl("dummy.cfg");
	BOOST_CHECK_EQUAL(data2.ab0c0d0->getIncludeUrl(), "dummy.cfg");
}

/**
 * Check if setRange and getRange work
 */
BOOST_AUTO_TEST_CASE(testRange)
{
	CConfig config(data2.ab0c0d0);
	config("x").setRange(Range::createFromString("[0.0, 10.0["));
	BOOST_CHECK_EQUAL(config("x").getRange()->toString().compare("[0.0, 10.0["), 0);
}

/**
 * Check if "at" method works (and indirectly appendPath).
 */
BOOST_AUTO_TEST_CASE(testAt)
{
	CConfig config(&data2.a);
	config.at("B[0].C[0].D[0].y");
	BOOST_CHECK_EQUAL(config.getRealValue(), -2.3);
}

///////////////////////////////////////////////////////////////////////////////
// fixture 3 test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check if search an object with 1 backtracking works.
 */
BOOST_AUTO_TEST_CASE(testSearch1Back3)
{
	CConfig config(&data3.cfg);
	Attribute* attr = config("A", "0")("B", "0")("C", "0")("D", "0").search();
	BOOST_CHECK_EQUAL(attr, data3.a0b_c_d0);
}

/**
 * Check if searching a scalar with nested backtrackings works.
 */
BOOST_AUTO_TEST_CASE(testSearchNestedBack3)
{
	CConfig config(&data3.cfg);
	Result res = config("A", "0")("B", "0")("C", "0")("D", "0")("x").resolve();
	BOOST_CHECK_EQUAL(res.getAttribute(), data3.a_b0c0d0x);
	BOOST_CHECK_EQUAL(res.getBacktrackLevel(), 3);
}

/**
 * Check if backtracking beyond the base object works.
 */
BOOST_AUTO_TEST_CASE(testSearchFullBack3)
{
	CConfig config(data3.a0b0c0);
	Attribute* attr = config("D", "0")("x").search();
	BOOST_CHECK_EQUAL(attr, data3.a_b0c0d0x);
}

/**
 * Check if using a config from default objects fails as expected.
 */
BOOST_AUTO_TEST_CASE(testSearchFail3)
{
	CConfig config(data3.a0b_c_d0);
	CATCH_ATTRIBUTE_NOT_FOUND("cfg.A[0].B[*].C[*].D[0].x",
		Attribute& attr = config("x").find();
		attr = 1
	);
}

/**
 * Check if config copy construction works.
 */
BOOST_AUTO_TEST_CASE(testCConfigCopyConstruct)
{
	CConfig config(data3.a0);

	CConfig config2 = config("B", "0")("C", "0")("D", "0");
	BOOST_CHECK_EQUAL(config2.getBaseAttribute(), data3.a0b_c_d0);

	Attribute* attr = config2("x").search();
	BOOST_CHECK_EQUAL(attr, data3.a_b0c0d0x);
}

/**
 * Check if config assignment works.
 */
BOOST_AUTO_TEST_CASE(testCConfigAssignment)
{
	CConfig config(data3.a0);
	CConfig config2(data3.a0b0);

	config2 = config("B", "0")("C", "0")("D", "0");
	BOOST_CHECK_EQUAL(config2.getBaseAttribute(), data3.a0b_c_d0);

	Attribute* attr = config2("x").search();
	BOOST_CHECK_EQUAL(attr, data3.a_b0c0d0x);
}

/**
 * Check if attribute enumeration occurs in the right order.
 */
BOOST_AUTO_TEST_CASE(testIterateOrder)
{
	CConfig config(&data3.cfg);
	CConfig::Iterator it = config("A", "0")("B", "0")("C", "0")("D", "0").iterate();

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: x");
	CConfig attr = it.next();
	BOOST_CHECK_EQUAL(attr.getBaseAttribute(), data3.a_b0c0d0x);
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 4);

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: y");
	it.next();
	attr = *it;
	BOOST_CHECK_EQUAL(attr.getBaseAttribute(), data3.a0b_c_d0y);
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 2);

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: z");
	attr = it.next();
	BOOST_CHECK_EQUAL(attr.getBaseAttribute(), data3.a_b0c0d0z);
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 4);

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: w");
	attr = it.next();
	BOOST_CHECK_EQUAL(attr.getBaseAttribute(), data3.a0b_c_d0w);
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 2);

	BOOST_CHECK_MESSAGE(not it.hasNext(), "no more attributes");
}

/**
 * Check if attribute enumeration occurs in the right order and if nested
 * enumeration works.
 */
BOOST_AUTO_TEST_CASE(testNestIterate)
{
	CConfig config(&data3.cfg);
	CConfig::Iterator it = config.iterate();

	// Browse cfg
	BOOST_CHECK_MESSAGE(it.hasNext(), "next: A[*]");
	CConfig attr = it.next();
	BOOST_CHECK_EQUAL(attr.getBaseAttribute(), data3.a_); // "A[*]"
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 0);

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: A[0]");
	CConfig a0 = it.next();
	BOOST_CHECK_EQUAL(a0.getBaseAttribute(), data3.a0); // A[0]
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 0);

	BOOST_CHECK_MESSAGE(not it.hasNext(), "no more attributes");

	// Browse A[0]
	it = a0.iterate();

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: B[0]");
	CConfig b0 = it.next();
	BOOST_CHECK_EQUAL(b0.getBaseAttribute(), data3.a0b0); // B[0]
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 0);

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: B[*]");
	attr = it.next();
	BOOST_CHECK_EQUAL(attr.getBaseAttribute(), data3.a0b_); // B[*]
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 0);

	BOOST_CHECK_MESSAGE(not it.hasNext(), "no more attributes");

	// Browse B[0]
	it = b0.iterate();

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: C[0]");
	CConfig c0 = it.next();
	BOOST_CHECK_EQUAL(c0.getBaseAttribute(), data3.a0b0c0); // C[0]
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 0);

	// C[*] should appear: inherited default objects are enumerated (since revision 1264).
	BOOST_CHECK_MESSAGE(it.hasNext(), "next: C[*]");
	attr = it.next();
	BOOST_CHECK_EQUAL(attr.getBaseAttribute(), data3.a0b_c_); // C[*]
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 1);

	BOOST_CHECK_MESSAGE(not it.hasNext(), "no more attributes");
}

/**
 * Check if scalar assignment of string works.
 */
BOOST_AUTO_TEST_CASE(testStringAssignment)
{
	CConfig config(&data3.cfg);
	config("A", "0")("B", "0")("C", "0")("D", "0")("x") = std::string("Hello");
	std::string x = data3.cfg("A", "0")("B", "0")("C", "0")("D", "0")("x");
	BOOST_CHECK_EQUAL(x, "Hello");
}

/**
 * Check if scalar assignment of int works.
 */
BOOST_AUTO_TEST_CASE(testIntAssignment)
{
	CConfig config(&data3.cfg);
	config("A", "0")("B", "0")("C", "0")("D", "0")("y") = 42;
	int y = data3.cfg("A", "0")("B", "0")("C", "0")("D", "0")("y");
	BOOST_CHECK_EQUAL(y, 42);
}

/**
 * Check if scalar assignment of hex works.
 */
BOOST_AUTO_TEST_CASE(testHexAssignment)
{
	CConfig config(&data3.cfg);
	config("A", "0")("B", "0")("C", "0")("D", "0")("z") = 0xBAADF00D;
	unsigned z = data3.cfg("A", "0")("B", "0")("C", "0")("D", "0")("z");
	BOOST_CHECK_EQUAL(z, 0xBAADF00D);
}

/**
 * Check if scalar assignment of real works.
 */
BOOST_AUTO_TEST_CASE(testRealAssignment)
{
	CConfig config(&data3.cfg);
	config("A", "0")("B", "0")("C", "0")("D", "0")("w") = -9.9;
	double w = data3.cfg("A", "0")("B", "0")("C", "0")("D", "0")("w");
	BOOST_CHECK_EQUAL(w, -9.9);
}

/**
 * Check if scalar assignment of boolean works.
 */
BOOST_AUTO_TEST_CASE(testBoolAssignment)
{
	CConfig config(&data3.cfg);
	config("A", "0")("B", "0")("C", "0")("D", "0")("x") = true;
	bool x = data3.cfg("A", "0")("B", "0")("C", "0")("D", "0")("x");
	BOOST_CHECK_EQUAL(x, true);
}

/**
 * Check if scalar assignment of const char* works.
 */
BOOST_AUTO_TEST_CASE(testCharAssignment)
{
	CConfig config(&data3.cfg);
	config("A", "0")("B", "0")("C", "0")("D", "0")("x") = "Hi There!";
	const char *x = data3.cfg("A", "0")("B", "0")("C", "0")("D", "0")("x");
	BOOST_CHECK(not strcmp(x, "Hi There!"));
}

/**
 * Check if string setter, getters and conversion operator works.
 */
BOOST_AUTO_TEST_CASE(testStringAccessors)
{
	std::string val = "Hello";
	Attribute& c = *data3.a0b0c0;
	CConfig config(&c);
	config("D", "0")("x").setStringValue(val);
	BOOST_CHECK_EQUAL(c("D", "0")("x").getStringValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getStringValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("absent").getStringValue(val), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getValueAsString(), val);
}

/**
 * Check if int setter, getters and conversion operator works.
 */
BOOST_AUTO_TEST_CASE(testIntAccessors)
{
	int val = 42;
	Attribute& c = *data3.a0b0c0;
	CConfig config(&c);

	config("D", "0")("x").setIntValue(val);
	BOOST_CHECK_EQUAL(c("D", "0")("x").getIntValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getIntValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("absent").getIntValue(val), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getValueAsInt(), val);

	config("D", "0")("x").setIntValue("4242");
	BOOST_CHECK_EQUAL(c("D", "0")("x").getIntValue(), 4242);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getIntValue(), 4242);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getValueAsInt(), 4242);
}

/**
 * Check if hex setter, getters and conversion operator works.
 */
BOOST_AUTO_TEST_CASE(testHexAccessors)
{
	unsigned val = 0xDEADBEEF;
	Attribute& c = *data3.a0b0c0;
	CConfig config(&c);

	config("D", "0")("x").setHexValue(val);
	BOOST_CHECK_EQUAL(c("D", "0")("x").getHexValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getHexValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("absent").getHexValue(val), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getValueAsHex(), val);

	config("D", "0")("x").setHexValue("0xABCD");
	BOOST_CHECK_EQUAL(c("D", "0")("x").getHexValue(), 0xABCDu);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getHexValue(), 0xABCDu);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getValueAsHex(), 0xABCDu);
}

/**
 * Check if real setter, getters and conversion operator works.
 */
BOOST_AUTO_TEST_CASE(testRealAccessors)
{
	double val = -9.9;
	Attribute& c = *data3.a0b0c0;
	CConfig config(&c);

	config("D", "0")("x").setRealValue(val);
	BOOST_CHECK_EQUAL(c("D", "0")("x").getRealValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getRealValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("absent").getRealValue(val), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getValueAsReal(), val);

	config("D", "0")("x").setRealValue("8.5");
	BOOST_CHECK_EQUAL(c("D", "0")("x").getRealValue(), 8.5);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getRealValue(), 8.5);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getValueAsReal(), 8.5);
}

/**
 * Check if bool setter, getters and conversion operator works.
 */
BOOST_AUTO_TEST_CASE(testBoolAccessors)
{
	bool val = true;
	Attribute& c = *data3.a0b0c0;
	CConfig config(&c);

	config("D", "0")("x").setBoolValue(val);
	BOOST_CHECK_EQUAL(c("D", "0")("x").getBoolValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getBoolValue(), val);
	BOOST_CHECK_EQUAL(config("D", "0")("absent").getBoolValue(val), val);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getValueAsBool(), val);

	config("D", "0")("x").setBoolValue(std::string("false"));
	BOOST_CHECK_EQUAL(c("D", "0")("x").getBoolValue(), false);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getBoolValue(), false);
	BOOST_CHECK_EQUAL(config("D", "0")("x").getValueAsBool(), false);
}

/**
 * Check if unit getter/setter works.
 */
BOOST_AUTO_TEST_CASE(testUnitAccessors)
{
	Attribute& c = *data3.a0b0c0;
	CConfig config(&c);
	config("D", "0")("x").setUnit("m");
	BOOST_CHECK_EQUAL(c("D", "0")("x").getUnit(), "m");
	BOOST_CHECK_EQUAL(config("D", "0")("x").getUnit(), "m");
}

/**
 * Check if getValueTypeId works.
 */
BOOST_AUTO_TEST_CASE(testGetValueTypeId)
{
	CConfig config(data3.a_b0c0d0x);
	BOOST_CHECK_EQUAL(config.getValueTypeId(), REAL);
}

/**
 * Check if create works.
 */
BOOST_AUTO_TEST_CASE(testCreate)
{
	CConfig config(data3.a0b0);
	Attribute& d = config("C", "0")("D", "0").create(OBJECT_INST);
	Attribute& d2 = config("C", "0")("D", "0").find();
	BOOST_CHECK_EQUAL(&d, &d2);
	BOOST_CHECK(typeid(d) == typeid(Object));
}

/**
 * Check create failure.
 */
BOOST_AUTO_TEST_CASE(testCreateFailre)
{
	/*try
	{
		CConfig config(data3.a_b0c0);
		Attribute& d = config("D", "0")("x")("y").create(SCALAR_INST);
	}*/
}

/**
 * Check if getAttributes works.
 */
BOOST_AUTO_TEST_CASE(testGetAttributes)
{
	CConfig config(&data3.cfg);
	ConfigList attributes = config("A", 0)("B", 0)("C", 0)("D", 0).getAttributes();
	BOOST_CHECK_EQUAL(attributes.size(), 4u);
	BOOST_CHECK_EQUAL(attributes[0].getPath().toString(), "cfg.A[0].B[0].C[0].D[0].x");
	BOOST_CHECK_EQUAL(attributes[1].getPath().toString(), "cfg.A[0].B[0].C[0].D[0].y");
	BOOST_CHECK_EQUAL(attributes[2].getPath().toString(), "cfg.A[0].B[0].C[0].D[0].z");
	BOOST_CHECK_EQUAL(attributes[3].getPath().toString(), "cfg.A[0].B[0].C[0].D[0].w");
}

/**
 * Check if setName works
 */
BOOST_AUTO_TEST_CASE(testCConfigSetName)
{
	CConfig config(data3.a0b0c0d_);
	config("x").setName("foo");
	BOOST_CHECK_EQUAL(config("foo").getName(), "foo");
	config.setName("E");
	BOOST_CHECK_EQUAL(config.getName(), "E");
	BOOST_CHECK_EQUAL(config.getPath().toString(), "cfg.A[0].B[0].C[0].E[*]");
}

/**
 * Check if setIndex works
 */
BOOST_AUTO_TEST_CASE(testSetIndex)
{
	CConfig config(data3.a0b0c0d_);
	config("x").setIndex("0");
	BOOST_CHECK_EQUAL(config("x", 0).getIndex(), "0");
	config.setIndex("bar");
	BOOST_CHECK_EQUAL(config.getIndex(), "bar");
	BOOST_CHECK_EQUAL(config.getPath().toString(), "cfg.A[0].B[0].C[0].D[bar]");
}

///////////////////////////////////////////////////////////////////////////////
// fixture 4 test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check if finding a default object without backtracking works.
 */
BOOST_AUTO_TEST_CASE(testSearchNoBack4)
{
	CConfig config(&data4.cfg);
	Attribute* attr = config("A", "*")("B", "*").search();
	BOOST_CHECK_EQUAL(attr, data4.a_b_);
}

/**
 * Check if finding a default object with backtracking works.
 */
BOOST_AUTO_TEST_CASE(testSearchBack4)
{
	CConfig config(data4.a0b0);
	Attribute* attr = config("C", "*").search();
	BOOST_CHECK_EQUAL(attr, data4.a0b_c_);
}

/**
 * Check if iteration loses information about the path when returning a
 * default object.
 */
BOOST_AUTO_TEST_CASE(testIteratePath)
{
	CConfig config(&data4.cfg);
	CConfig::Iterator it = config("A", 0)("B", 0).iterate();

	// Browse A[0].B[0]
	// C[*] should appear: inherited default objects are enumerated (since revision 1264).

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: C[*]");
	CConfig attr = it.next();
	BOOST_CHECK_EQUAL(attr.getBaseAttribute(), data4.a0b_c_); // C[*]
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 1);

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: C[0]");
	CConfig c0 = it.next();
	BOOST_CHECK_EQUAL(c0.getBaseAttribute(), data4.a0b_c0); // C[0]
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 1);

	BOOST_CHECK_MESSAGE(not it.hasNext(), "no more attributes");

	// Browse virtual object A[0].B[0].C[0]
	it = c0.iterate();

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: u");
	Attribute* u = it.next().search();
	BOOST_CHECK_EQUAL(u, data4.a_b0c_u); // C[*]
	BOOST_CHECK_EQUAL(it.get().getBacktrackLevel(), 3);

	BOOST_CHECK_MESSAGE(not it.hasNext(), "no more attributes");
}

///////////////////////////////////////////////////////////////////////////////
// fixture 5 test methods body
///////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(testIterate5)
{
	CConfig config(&data5.cfg);
	CConfig::Iterator it = config("A", 0).iterate("c");

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: c[0]");
	CConfig c0 = it.next();
	BOOST_CHECK_EQUAL(c0.getBaseAttribute(), data5.a_c0); // c[0]

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: c[1]");
	CConfig c1 = it.next();
	BOOST_CHECK_EQUAL(c1.getBaseAttribute(), data5.a0c1); // c[1]

	BOOST_CHECK_MESSAGE(it.hasNext(), "next: c[2]");
	CConfig c2 = it.next();
	BOOST_CHECK_EQUAL(c2.getBaseAttribute(), data5.a0c2); // c[2]

	BOOST_CHECK_MESSAGE(not it.hasNext(), "no more attributes");
}

BOOST_AUTO_TEST_SUITE_END()

}
