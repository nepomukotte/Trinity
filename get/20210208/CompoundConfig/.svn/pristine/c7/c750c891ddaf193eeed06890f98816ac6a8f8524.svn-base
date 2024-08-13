// =============================================================================
// InheritersIteratorTest.cpp
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

#include "CCfg/View/InheritersIterator.h"
#include "CCfg/CConfig.h"
#include "CCfg/Object.h"
#include "CCfg/Attribute.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{
namespace View
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_view_heir
{
	Object cfg;

	Object *a_;
	Object *a_b_;
	Object *a_b0;
	Object *a_b1;

	Object *a0;
	Object *a0b_;
	Object *a0b0;
	Object *a0b2;

	Object *a1;
	Object *a1b3;

	config_view_heir()
		: cfg("cfg")
	{
		a_ = new Object("A", "*");
		a_b_ = new Object("B", "*");
		a_b0 = new Object("B", "0");
		a_b1 = new Object("B", "1");

		a0 = new Object("A", "0");
		a0b_ = new Object("B", "*");
		a0b0 = new Object("B", "0");
		a0b2 = new Object("B", "2");

		a1 = new Object("A", "1");
		a1b3 = new Object("B", "3");

		cfg.append(a_);
		a_->append(a_b_);
		a_->append(a_b0);
		a_->append(a_b1);
		cfg.append(a0);
		a0->append(a0b_);
		a0->append(a0b0);
		a0->append(a0b2);
		cfg.append(a1);
		a1->append(a1b3);
	}

	~config_view_heir()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_View_InheritersIterator, config_view_heir)

///////////////////////////////////////////////////////////////////////////////
// fixture 1 test methods body
///////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(testInheritersIterator1)
{
	CConfig config(&cfg);
	InheritersIterator heir = config("A", "*")("B", "*").getInheriterAttributes();

	BOOST_CHECK_MESSAGE(heir.hasNext(), "next1");
	heir.next();
	BOOST_CHECK_EQUAL(heir.get().getAttribute(), a_b0);

	BOOST_CHECK_MESSAGE(heir.hasNext(), "next2");
	heir.next();
	BOOST_CHECK_EQUAL(heir.get().getAttribute(), a_b1);

	BOOST_CHECK_MESSAGE(heir.hasNext(), "next3");
	heir.next();
	BOOST_CHECK_EQUAL(heir.get().getAttribute(), a0b_);

	BOOST_CHECK_MESSAGE(heir.hasNext(), "next4");
	heir.next();
	BOOST_CHECK_EQUAL(heir.get().getAttribute(), a0b0);

	BOOST_CHECK_MESSAGE(heir.hasNext(), "next5");
	heir.next();
	BOOST_CHECK_EQUAL(heir.get().getAttribute(), a0b2);

	BOOST_CHECK_MESSAGE(heir.hasNext(), "next6");
	heir.next();
	BOOST_CHECK_EQUAL(heir.get().getAttribute(), a1b3);

	BOOST_CHECK_MESSAGE(not heir.hasNext(), "no more objects");
}

BOOST_AUTO_TEST_CASE(testInheritersIterator2)
{
	CConfig config(&cfg);
	InheritersIterator heir = config("A", "*")("B", "0").getInheriterAttributes();

	BOOST_CHECK_MESSAGE(heir.hasNext(), "next1");
	heir.next();
	BOOST_CHECK_EQUAL(heir.get().getAttribute(), a0b0);

	BOOST_CHECK_MESSAGE(not heir.hasNext(), "no more objects");
}

BOOST_AUTO_TEST_CASE(testInheritersIterator3)
{
	CConfig config(&cfg);
	InheritersIterator heir = config("A", "0")("B", "2").getInheriterAttributes();

	BOOST_CHECK_MESSAGE(not heir.hasNext(), "no more objects");
}

BOOST_AUTO_TEST_SUITE_END()

}
}
