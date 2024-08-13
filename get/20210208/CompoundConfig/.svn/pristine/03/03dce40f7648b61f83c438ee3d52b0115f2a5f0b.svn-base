// =============================================================================
// IntervalTest.cpp
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

#include "CCfg/Interval.h"
#include <boost/test/auto_unit_test.hpp>
#include <limits>

#undef max

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_interval
{
	config_interval()
	{
	}

	~config_interval()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Interval, config_interval)


///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Test closed integer intervals.
 */
BOOST_AUTO_TEST_CASE(testClosedIntegerIntervals)
{
	Interval interval(-1, 1);
	BOOST_CHECK(interval.contains(Value(-1)));
	BOOST_CHECK(interval.contains(Value(0)));
	BOOST_CHECK(interval.contains(Value(1)));
	BOOST_CHECK(not interval.contains(Value(-2)));
	BOOST_CHECK(not interval.contains(Value(2)));
	BOOST_CHECK(not interval.contains(Value(0.5)));

	std::string str = interval.toString();
	BOOST_CHECK_EQUAL(str, "[-1, 1]");

	Interval parsedInterval;
	parsedInterval.fromString(str);
	BOOST_CHECK_EQUAL(parsedInterval, interval);
}

/**
 * Test open integer intervals.
 */
BOOST_AUTO_TEST_CASE(testOpenedIntegerIntervals)
{
	Interval interval(-1, 1, Interval::ExcludeBound, Interval::ExcludeBound);
	BOOST_CHECK(not interval.contains(Value(-1)));
	BOOST_CHECK(interval.contains(Value(0)));
	BOOST_CHECK(not interval.contains(Value(1)));
	BOOST_CHECK(not interval.contains(Value(-2)));
	BOOST_CHECK(not interval.contains(Value(2)));
	BOOST_CHECK(not interval.contains(Value(-0.5)));

	std::string str = interval.toString();
	BOOST_CHECK_EQUAL(str, "]-1, 1[");

	Interval parsedInterval;
	parsedInterval.fromString(str);
	BOOST_CHECK_EQUAL(parsedInterval, interval);
}

/**
 * Test hexadecimal intervals.
 */
BOOST_AUTO_TEST_CASE(testHexIntervals)
{
	Interval interval(0xF0u, 0x100u, Interval::IncludeBound, Interval::ExcludeBound);
	BOOST_CHECK(interval.contains(Value(0xF0u)));
	BOOST_CHECK(interval.contains(Value(0xFAu)));
	BOOST_CHECK(interval.contains(Value(0xFFu)));
	BOOST_CHECK(not interval.contains(Value(0xE0u)));
	BOOST_CHECK(not interval.contains(Value(0x100u)));
	BOOST_CHECK(not interval.contains(Value(256)));

	std::string str = interval.toString();
	BOOST_CHECK_EQUAL(str, "[0xF0, 0x100[");

	Interval parsedInterval;
	parsedInterval.fromString(str);
	BOOST_CHECK_EQUAL(parsedInterval, interval);
}

/**
 * Test floatting point intervals.
 */
BOOST_AUTO_TEST_CASE(testRealIntervals)
{
	Interval interval(-3.3, 3.3, Interval::IncludeBound, Interval::ExcludeBound);
	BOOST_CHECK(interval.contains(Value(-3.3)));
	BOOST_CHECK(interval.contains(Value(1.5)));
	BOOST_CHECK(not interval.contains(Value(3.3)));
	BOOST_CHECK(not interval.contains(Value(3)));

	// *2 is for floor(log2(3.3)), because epsilon has to be put on the same scale as the number.
	BOOST_CHECK(interval.contains(Value(3.3 - std::numeric_limits<double>::epsilon()*2)));
	BOOST_CHECK(not interval.contains(Value(3.3 - std::numeric_limits<double>::epsilon())));

	std::string str = interval.toString();
	BOOST_CHECK_EQUAL(str, "[-3.3, 3.3[");

	Interval parsedInterval;
	parsedInterval.fromString(str);
	BOOST_CHECK_EQUAL(parsedInterval, interval);
}

/**
 * Test string intervals.
 */
BOOST_AUTO_TEST_CASE(testStringIntervals)
{
	Interval interval("ar", "as", Interval::ExcludeBound, Interval::IncludeBound);
	BOOST_CHECK(not interval.contains(Value("ar")));
	BOOST_CHECK(interval.contains(Value("art")));
	BOOST_CHECK(interval.contains(Value("as")));
	BOOST_CHECK(not interval.contains(Value("aquarium")));
	BOOST_CHECK(not interval.contains(Value("assistant")));
	BOOST_CHECK(not interval.contains(Value(42)));

	std::string str = interval.toString();
	BOOST_CHECK_EQUAL(str, "]ar, as]");

	Interval parsedInterval;
	parsedInterval.fromString(str);
	BOOST_CHECK_EQUAL(parsedInterval, interval);
}

/**
 * Test escaped string intervals.
 */
BOOST_AUTO_TEST_CASE(testEscStringIntervals)
{
	Interval interval("A[0].B[0]", "A[0].B[0^", Interval::IncludeBound, Interval::ExcludeBound);
	BOOST_CHECK(interval.contains(Value("A[0].B[0]")));
	BOOST_CHECK(interval.contains(Value("A[0].B[0].x")));
	BOOST_CHECK(not interval.contains(Value("A[0].A[0]")));
	BOOST_CHECK(not interval.contains(Value("A[0].B[1]")));
	BOOST_CHECK(not interval.contains(Value(42)));

	std::string str = interval.toString();
	BOOST_CHECK_EQUAL(str, "[\"A[0].B[0]\", \"A[0].B[0^\"[");

	Interval parsedInterval;
	parsedInterval.fromString(str);
	BOOST_CHECK_EQUAL(parsedInterval, interval);
}

/**
 * Test equality comparisons
 */
BOOST_AUTO_TEST_CASE(testEquals)
{
	Interval interval(-2.5, 2.5, Interval::IncludeBound, Interval::ExcludeBound);
	BOOST_CHECK(interval.equals(interval));
	BOOST_CHECK(interval == interval);

	Interval interval2(-2.5, 2.5);
	BOOST_CHECK(not interval.equals(interval2));
	BOOST_CHECK(interval != interval2);

	Interval interval3(-2.5, 3.0);
	BOOST_CHECK(not interval.equals(interval3));
	BOOST_CHECK(interval != interval3);

	Interval interval4(-2, 2, Interval::IncludeBound, Interval::ExcludeBound);
	BOOST_CHECK(not interval.equals(interval4));
	BOOST_CHECK(interval != interval4);
}

/**
 * Test infinity interval on both bounds
 */
BOOST_AUTO_TEST_CASE(testInfinityIntervals)
{
	Interval interval(0.0, 0.0, Interval::IncludeBound, Interval::Unbounded);
	BOOST_CHECK(interval.contains(Value(0.0)));
	BOOST_CHECK(interval.contains(Value(10.0)));
	BOOST_CHECK(interval.contains(Value(std::numeric_limits<double>::max())));
	BOOST_CHECK(not interval.contains(Value(-1.0)));
	BOOST_CHECK(not interval.contains(Value("hello")));

	std::string str = interval.toString();
	BOOST_CHECK_EQUAL(str, "[0.0, +inf[");

	Interval parsedInterval;
	parsedInterval.fromString(str);
	BOOST_CHECK_EQUAL(parsedInterval, interval);
}

/**
 * Test infinity interval on both bounds
 */
BOOST_AUTO_TEST_CASE(testInfinityIntervals2)
{
	Interval interval(0, 0, Interval::Unbounded, Interval::Unbounded);
	BOOST_CHECK(interval.contains(Value(-std::numeric_limits<double>::max())));
	BOOST_CHECK(interval.contains(Value(std::numeric_limits<double>::max())));
	BOOST_CHECK(interval.contains(Value(std::numeric_limits<int>::max())));
	BOOST_CHECK(interval.contains(Value(std::numeric_limits<unsigned>::max())));
	BOOST_CHECK(interval.contains(Value("hello")));

	std::string str = interval.toString();
	BOOST_CHECK_EQUAL(str, "]-inf, +inf[");

	Interval parsedInterval;
	parsedInterval.fromString(str);
	BOOST_CHECK_EQUAL(parsedInterval.toString(), interval.toString());
	BOOST_CHECK_EQUAL(parsedInterval, interval);
}

BOOST_AUTO_TEST_SUITE_END()

}
