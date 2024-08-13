// =============================================================================
// RangeTest.cpp
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

#include "CCfg/Range.h"
#include "CCfg/Interval.h"
#include "CCfg/Enum.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_range
{
	config_range()
	{
	}

	~config_range()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Range, config_range)


///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Test enumeration parsing.
 */
BOOST_AUTO_TEST_CASE(testEnumParsing)
{
	RangePtr r = Range::createFromString("{-2.5, 0.0, 1.0}");
	BOOST_CHECK(dynamic_cast<Enum*>(r.get()));
	BOOST_CHECK_EQUAL(r->toString(), "{-2.5, 0.0, 1.0}");
}

/**
 * Test interval parsing.
 */
BOOST_AUTO_TEST_CASE(testIntervalParsing)
{
	RangePtr r = Range::createFromString("[-1.0, 1.0[");
	BOOST_CHECK(dynamic_cast<Interval*>(r.get()));
	BOOST_CHECK_EQUAL(r->toString(), "[-1.0, 1.0[");
}

/**
 * Test invalid range parsing.
 */
BOOST_AUTO_TEST_CASE(testRangeParsing)
{
	RangePtr r = Range::createFromString("[-1.0, 1.0, 5.0[");
	BOOST_CHECK(r.get() == 0);
}

/**
 * Test invalid range parsing.
 */
BOOST_AUTO_TEST_CASE(testRangeParsingInvalid)
{
	RangePtr r = Range::createFromString("hello");
	BOOST_CHECK(r.get() == 0);
}

BOOST_AUTO_TEST_SUITE_END()

}
