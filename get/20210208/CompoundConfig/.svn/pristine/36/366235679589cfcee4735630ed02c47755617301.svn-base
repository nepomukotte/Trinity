// =============================================================================
// EnumTest.cpp
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

#include "CCfg/Enum.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_enum
{
	config_enum()
	{
	}

	~config_enum()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Enum, config_enum)


///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Test initialization, concatenation, reinitialization for each type.
 */
BOOST_AUTO_TEST_CASE(testInitConcatReinit)
{
	Enum e;
	e = 1, 8, 3;
	BOOST_CHECK_EQUAL(e.toString(), "{1, 3, 8}");
	e += 12, 10;
	BOOST_CHECK_EQUAL(e.toString(), "{1, 3, 8, 10, 12}");

	e = 0xAu, 0xFu, 0xFFu;
	BOOST_CHECK_EQUAL(e.toString(), "{0xA, 0xF, 0xFF}");

	e = 3.3, 5.0, 12.0;
	BOOST_CHECK_EQUAL(e.toString(), "{3.3, 5.0, 12.0}");

	e = true, false;
	BOOST_CHECK_EQUAL(e.toString(), "{false, true}");

	e = "hello", "world", "1,000,000";
	BOOST_CHECK_EQUAL(e.toString(), "{\"1,000,000\", hello, world}");

	e = "\x7", "\t", "\"hi!\"";
	BOOST_CHECK_EQUAL(e.toString(), "{\"\\x7\", \"\\t\", \"\\\"hi!\\\"\"}");
}

/**
 * Test contains.
 */
BOOST_AUTO_TEST_CASE(testContains)
{
	Enum e;
	e = 1, 3, 8;
	BOOST_CHECK(e.contains(Value(1)));
	BOOST_CHECK(e.contains(Value(8)));
	BOOST_CHECK(not e.contains(Value(2)));

	e = "hello", "world";
	BOOST_CHECK(e.contains(Value("hello")));
	BOOST_CHECK(e.contains(Value("world")));
	BOOST_CHECK(not e.contains(Value("blue")));
}

/**
 * Test copy constructor and assignment.
 */
BOOST_AUTO_TEST_CASE(testCConstructorAssignment)
{
	Enum e;
	e = 1, 5, 10;
	Enum e1 = e;
	BOOST_CHECK_EQUAL(e.toString(), "{1, 5, 10}");

	e = 10, 20, 30;
	e1 = e;
	BOOST_CHECK_EQUAL(e.toString(), "{10, 20, 30}");
}

/**
 * Test linear list
 */
BOOST_AUTO_TEST_CASE(testLinearList)
{
	Enum e;
	e = 0, Enum::linear(100, 110, 5), 1000;
	BOOST_CHECK_EQUAL(e.toString(), "{0, 100, 105, 110, 1000}");

	e = Enum::linear(0x0u, 0x10u, 4u), 0xFF;
	BOOST_CHECK_EQUAL(e.toString(), "{0x0, 0x4, 0x8, 0xC, 0x10, 0xFF}");

	e += Enum::linear(0x20u, 0x30u, 0x10u);
	BOOST_CHECK_EQUAL(e.toString(), "{0x0, 0x4, 0x8, 0xC, 0x10, 0x20, 0x30, 0xFF}");
}

/**
 * Test equality comparisons
 */
BOOST_AUTO_TEST_CASE(testEquality)
{
	Enum e1;
	e1 = 0.0, 1.0/4, 1.0/2, 3.0/4, 1.0;
	Enum e2;
	e2 = Enum::linear(0.0, 1.0, 1.0/4);
	BOOST_CHECK_EQUAL(e1, e2);
	BOOST_CHECK(e1.equals(e2));

	e2.insert(5.0/4);
	BOOST_CHECK(e1 != e2);
}

/**
 * Test type homogeneity exceptions
 */
BOOST_AUTO_TEST_CASE(testTypeError)
{
	try
	{
		Enum e;
		e = 1, 2;
		e += 3.5;
		BOOST_FAIL("Should throw IncompatibleType");
	}
	catch(Enum::IncompatibleType& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"There was an attempt to insert the value: 3.5 (real number) "
			"in the enum {1, 2} of type \"integer\".");
	}
}

/**
 * Test fromString
 */
BOOST_AUTO_TEST_CASE(testEnumFromString)
{
	Enum e;
	e = 1;
	e.fromString("{}");
	BOOST_CHECK_EQUAL(e.toString(), "{}");

	e.fromString("{ 50 }");
	BOOST_CHECK_EQUAL(e.toString(), "{50}");
	BOOST_CHECK_EQUAL(e.getValues().begin()->getType(), INTEGER);

	e.fromString("{ 3 , 1, 8  }");
	BOOST_CHECK_EQUAL(e.toString(), "{1, 3, 8}");
	BOOST_CHECK_EQUAL(e.getValues().begin()->getType(), INTEGER);

	e.fromString("{0xA, 0xFF, 0xF}");
	BOOST_CHECK_EQUAL(e.toString(), "{0xA, 0xF, 0xFF}");
	BOOST_CHECK_EQUAL(e.getValues().begin()->getType(), HEXADECIMAL);

	e.fromString("{3.3, 5.0, 12.0}");
	BOOST_CHECK_EQUAL(e.toString(), "{3.3, 5.0, 12.0}");
	BOOST_CHECK_EQUAL(e.getValues().begin()->getType(), REAL);

	e.fromString("{true, false}");
	BOOST_CHECK_EQUAL(e.toString(), "{false, true}");
	BOOST_CHECK_EQUAL(e.getValues().begin()->getType(), BOOLEAN);

	e.fromString("{  \"1,000,000\" ,  hello  ,  \"hello world !!\" }");
	BOOST_CHECK_EQUAL(e.toString(), "{\"1,000,000\", hello, \"hello world !!\"}");
	BOOST_CHECK_EQUAL(e.getValues().begin()->getType(), STRING);

	e.fromString("{\"\\x7\", \"\\t\", \"\\\"hi!\\\"\"}");
	BOOST_CHECK_EQUAL(e.toString(), "{\"\\x7\", \"\\t\", \"\\\"hi!\\\"\"}");

	try
	{
		e.fromString("{1, 5.0}");
		BOOST_FAIL("Should throw IncompatibleType exception");
	}
	catch(Enum::IncompatibleType& e)
	{
	}
}

BOOST_AUTO_TEST_SUITE_END()

}
