// =============================================================================
// ClassIteratorTest.cpp
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
#include "CCfg/View/ClassIterator.h"
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

struct config_view_classiterator
{
	Config3 data3;

	config_view_classiterator()
	{

	}

	~config_view_classiterator()
	{

	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_View_ClassIterator, config_view_classiterator)

///////////////////////////////////////////////////////////////////////////////
// fixture 3 test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 *
 */
BOOST_AUTO_TEST_CASE(testClassIterator1)
{
	CConfig config(&data3.cfg);
	ClassIterator iInst = config("A", 0)("B", 0).getClassInstances();

	BOOST_CHECK_MESSAGE(iInst.hasNext(), "next1");
	iInst.next();
	BOOST_CHECK_EQUAL(iInst.get().getAttribute(), data3.a_b0);

	BOOST_CHECK_MESSAGE(iInst.hasNext(), "next2");
	iInst.next();
	BOOST_CHECK_EQUAL(iInst.get().getAttribute(), data3.a0b_);

	BOOST_CHECK_MESSAGE(iInst.hasNext(), "next3");
	iInst.next();
	BOOST_CHECK_EQUAL(iInst.get().getAttribute(), data3.a0b0);

	BOOST_CHECK_MESSAGE(not iInst.hasNext(), "no more objects");
}

/**
 *
 */
BOOST_AUTO_TEST_CASE(testClassIterator2)
{
	CConfig config(&data3.cfg);
	ClassIterator iInst = config("A", 0)("B", 0)("C", 0)("D", 0)("y").getClassInstances();

	BOOST_CHECK_MESSAGE(iInst.hasNext(), "next1");
	iInst.next();
	BOOST_CHECK_EQUAL(iInst.get().getAttribute(), data3.a_b0c0d0y);

	BOOST_CHECK_MESSAGE(iInst.hasNext(), "next2");
	iInst.next();
	BOOST_CHECK_EQUAL(iInst.get().getAttribute(), data3.a0b_c_d0y);

	BOOST_CHECK_MESSAGE(not iInst.hasNext(), "no more attributes");
}

/**
 *
 */
BOOST_AUTO_TEST_CASE(testClassIterator3)
{
	CConfig config(&data3.cfg);
	ClassIterator iInst = config("A", 0)("B", 0)("C", 0)("D", 0)("u").getClassInstances();

	BOOST_CHECK_MESSAGE(not iInst.hasNext(), "no more attributes");
}

/**
 *
 */
BOOST_AUTO_TEST_CASE(testClassIterator4)
{
	ClassIterator iInst(data3.cfg, Path("cfg.A[0].B[0].C[0].D[0].x"));

	BOOST_CHECK_MESSAGE(iInst.hasNext(), "next1");
	iInst.next();
	BOOST_CHECK_EQUAL(iInst.get().getAttribute(), data3.a_b0c0d0x);

	BOOST_CHECK_MESSAGE(iInst.hasNext(), "next2");
	iInst.next();
	BOOST_CHECK_EQUAL(iInst.get().getAttribute(), data3.a0b0c0d_x);

	BOOST_CHECK_MESSAGE(not iInst.hasNext(), "no more attributes");
}

/**
 *
 */
BOOST_AUTO_TEST_CASE(testClassIterator5)
{
	ClassIterator iInst(*data3.a0b_c_d0w);

	BOOST_CHECK_MESSAGE(iInst.hasNext(), "next1");
	iInst.next();
	BOOST_CHECK_EQUAL(iInst.get().getAttribute(), data3.a0b_c_d0w);

	BOOST_CHECK_MESSAGE(not iInst.hasNext(), "no more attributes");
}

BOOST_AUTO_TEST_SUITE_END()

}
}
