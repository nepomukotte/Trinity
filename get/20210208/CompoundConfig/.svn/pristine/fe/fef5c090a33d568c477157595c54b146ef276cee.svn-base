// =============================================================================
// InheritedIteratorTest.cpp
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

#include "Configs.h"
#include "CCfg/View/InheritedIterator.h"
#include "CCfg/CConfig.h"
#include "CCfg/Path.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{
namespace View
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_view_default
{
	Config1 data1;
	Config2 data2;
	Config3 data3;

	config_view_default()
	{

	}

	~config_view_default()
	{

	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_View_InheritedIterator, config_view_default)

///////////////////////////////////////////////////////////////////////////////
// fixture 1 test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check if enumerating the default objects of A[0] works.
 */
BOOST_AUTO_TEST_CASE(testInheritedIterator1)
{
	InheritedIterator def(*data1.a1);
	BOOST_CHECK_MESSAGE(def.hasNext(), "next1");
	def.next();
	BOOST_CHECK_EQUAL(def.get().getAttribute(), data1.a_);
	BOOST_CHECK_EQUAL(def.get().getBacktrackLevel(), 1);

	BOOST_CHECK_MESSAGE(not def.hasNext(), "no more objects");
}

///////////////////////////////////////////////////////////////////////////////
// fixture 2 test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check if enumerating the default objects of an object that doesn't have any
 * works.
 */
BOOST_AUTO_TEST_CASE(testInheritedIterator2)
{
	InheritedIterator def(data2.a);

	BOOST_CHECK_MESSAGE(not def.hasNext(), "no objects");
}

/**
 * Check if enumerating the default objects with backtracking works.
 */
BOOST_AUTO_TEST_CASE(testInheritedIterator3)
{
	InheritedIterator def(*data2.ab0c0);

	BOOST_CHECK_MESSAGE(def.hasNext(), "next1");
	def.next();
	BOOST_CHECK_EQUAL(def.get().getAttribute(), data2.ab_c_);
	BOOST_CHECK_EQUAL(def.get().getBacktrackLevel(), 2);

	BOOST_CHECK_MESSAGE(not def.hasNext(), "no more objects");
}

/**
 * Check if enumerating default objects at different backtracking levels works.
 */
BOOST_AUTO_TEST_CASE(testInheritedIterator4)
{
	InheritedIterator def(*data2.ab0c0d0);

	BOOST_CHECK_MESSAGE(def.hasNext(), "next1");
	def.next();
	BOOST_CHECK_EQUAL(def.get().getAttribute(), data2.ab0c0d_);
	BOOST_CHECK_EQUAL(def.get().getBacktrackLevel(), 1);

	BOOST_CHECK_MESSAGE(def.hasNext(), "next2");
	def.next();
	BOOST_CHECK_EQUAL(def.get().getAttribute(), data2.ab_c_d0);
	BOOST_CHECK_EQUAL(def.get().getBacktrackLevel(), 3);

	BOOST_CHECK_MESSAGE(not def.hasNext(), "no more objects");
}

/**
 * Check if enumerating default objects of a default object fails as expected.
 */
BOOST_AUTO_TEST_CASE(testInheritedIterator5)
{
	InheritedIterator def(*data2.ab0c0d_);

	BOOST_CHECK_MESSAGE(not def.hasNext(), "no objects");
}

///////////////////////////////////////////////////////////////////////////////
// fixture 3 test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check if enumerating default values of a scalar works.
 */
BOOST_AUTO_TEST_CASE(testInheritedIterator6)
{
	CConfig config(&data3.cfg);
	InheritedIterator def = config("A", 0)("B", 0)("C", 0)("D", 0)("y").getInheritedAttributes();

	BOOST_CHECK_MESSAGE(def.hasNext(), "next1");
	def.next();
	BOOST_CHECK_EQUAL(def.get().getAttribute(), data3.a_b0c0d0y);
	BOOST_CHECK_EQUAL(def.get().getBacktrackLevel(), 5);

	BOOST_CHECK_MESSAGE(not def.hasNext(), "no more objects");
}

BOOST_AUTO_TEST_CASE(testInheritedIterator7)
{
	CConfig config(&data3.cfg);
	InheritedIterator def = config("A", 0)("B", 0)("C", 0)("D", 0).getInheritedAttributes();

	BOOST_CHECK_MESSAGE(def.hasNext(), "next1");
	def.next();
	BOOST_CHECK_EQUAL(def.get().getAttribute(), data3.a_b0c0d0);
	BOOST_CHECK_EQUAL(def.get().getBacktrackLevel(), 4);

	BOOST_CHECK_MESSAGE(not def.hasNext(), "no more objects");
}

BOOST_AUTO_TEST_SUITE_END()

}
}
