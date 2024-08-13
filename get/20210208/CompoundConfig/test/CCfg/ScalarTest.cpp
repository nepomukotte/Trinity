// =============================================================================
// ScalarTest.cpp
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

#include "CCfg/Object.h"
#include "CCfg/Scalar.h"
#include "CCfg/ClassInfo.h"
#include "CCfg/Enum.h"
#include "CCfg/Interval.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_scalar_data
{
	Object config;
	Object* a_;
	Object* a0;
	Scalar* a_x;
	Scalar* a_y;
	Scalar* a0x;
	Scalar* a0x1;

	config_scalar_data()
		: config("Daq")
	{
		a_ = new Object("A", "*");
		a0 = new Object("A", "0");
		a_x = new Scalar("x");
		a_y = new Scalar("y");
		a0x = new Scalar("x");
		a0x1 = new Scalar("x", "1");
		config.append(a_);
		config.append(a0);
		a_->append(a_x);
		a_->append(a_y);
		a0->append(a0x);
		a0->append(a0x1);

		a_y->setUnit("m");
	}

	~config_scalar_data()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Scalar, config_scalar_data)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check unit propagation
 */
BOOST_AUTO_TEST_CASE(testUnitPropagation)
{
	a_x->setUnit("V");
	BOOST_CHECK_EQUAL(a0x->getUnit(), "V");
	BOOST_CHECK_EQUAL(a0x1->getUnit(), "V");
	BOOST_CHECK_EQUAL(a_y->getUnit(), "m");
	BOOST_CHECK_EQUAL(a0x1->getClassInfo(), a0x->getClassInfo());
	BOOST_CHECK_EQUAL(a0x->getClassInfo(), a_x->getClassInfo());

	a0x->setUnit("mV");
	BOOST_CHECK_EQUAL(a_x->getUnit(), "mV");
}

/**
 * Check range propagation
 */
BOOST_AUTO_TEST_CASE(testRangePropagation)
{
	Enum range;
	range = 1, 2, 3;
	a0x1->setRange(range);

	BOOST_CHECK(*a_x->getRange() == range);
}

/**
 * Check ClassInfo attachment after setName.
 */
BOOST_AUTO_TEST_CASE(testScalarSetName)
{
	a0x1->setName("y");
	BOOST_CHECK_EQUAL(a0x1->getClassInfo(), a_y->getClassInfo());
	BOOST_CHECK_EQUAL(a0x->getClassInfo(), a_x->getClassInfo());
}

/**
 * Check ClassInfo attachment failure after setName (case destination has no classinfo).
 */
BOOST_AUTO_TEST_CASE(testSetNameFailure1)
{
	a_y->setIndex("0"); // needed to avoid name conflict.
	try
	{
		a_y->setName("x");
	}
	catch(Attribute::IncompatibleClassInfo& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The informations about the configuration class of the attribute "
			"Daq.A[*].x[0] which are: {unit=\"m\" range=null} "
			"are incompatible with those of the class it has joined: null");
	}
}

/**
 * Check ClassInfo attachment failure after setName (case destination unit different).
 */
BOOST_AUTO_TEST_CASE(testSetNameFailure2)
{
	a_x->setUnit("km");
	a_y->setIndex("0"); // needed to avoid name conflict.
	try
	{
		a_y->setName("x");
	}
	catch(Attribute::IncompatibleClassInfo& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The informations about the configuration class of the attribute "
			"Daq.A[*].x[0] which are: {unit=\"m\" range=null} "
			"are incompatible with those of the class it has joined: "
			"{unit=\"km\" range=null}");
	}
}

/**
 * Check ClassInfo attachment failure after setName (case destination range different).
 */
BOOST_AUTO_TEST_CASE(testSetNameFailure3)
{
	a_x->setUnit("m");
	Interval i(0, 10);
	a_x->setRange(i);
	a_y->setIndex("0"); // needed to avoid name conflict.
	try
	{
		a_y->setName("x");
	}
	catch(Attribute::IncompatibleClassInfo& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"The informations about the configuration class of the attribute "
			"Daq.A[*].x[0] which are: {unit=\"m\" range=null} "
			"are incompatible with those of the class it has joined: "
			"{unit=\"m\" range=[0, 10]}");
	}
}

/**
 * Check ClassInfo attachment success after setName when unit are the same.
 */
BOOST_AUTO_TEST_CASE(testClassInfoAttachment)
{
	a_x->setUnit("m"); // we will need a ClassInfo for x
	a_y->setIndex("0"); // needed to avoid name conflict.

	// preconditions
	BOOST_CHECK(a_x->getClassInfo() != a_y->getClassInfo());
	BOOST_CHECK_EQUAL(*a_x->getClassInfo(), *a_y->getClassInfo());

	// rename
	a_y->setName("x");

	// postconditions
	BOOST_CHECK_EQUAL(a_x->getClassInfo(), a_y->getClassInfo()); // class info have merged
}

/**
 * Check ClassInfo detachment when creating a new class with setName.
 */
BOOST_AUTO_TEST_CASE(testClassInfoDetachment)
{
	a0x1->setUnit("m"); // we will need a ClassInfo for x

	// preconditions: ensure the same ClassInfo is shared
	BOOST_CHECK_EQUAL(a0x1->getClassInfo(), a_x->getClassInfo());

	// change the name for an unexisting name
	a0x1->setName("z");

	// postconditions: ensure the ClassInfo has been cloned.
	BOOST_CHECK(a0x1->getClassInfo() != a_x->getClassInfo());
	BOOST_CHECK_EQUAL(*a0x1->getClassInfo(), *a_x->getClassInfo());
}

/**
 * Check ClassInfo detachment when creating a new class with moveTo.
 */
BOOST_AUTO_TEST_CASE(testScalarMoveTo)
{
	a0x1->setUnit("m"); // we will need a ClassInfo for x
	Object* b0 = new Object("B", "0");
	config.append(b0);

	// preconditions: ensure the same ClassInfo is shared
	BOOST_CHECK_EQUAL(a0x1->getClassInfo(), a_x->getClassInfo());

	// move x[1]
	a0x1->moveTo(*b0);

	// postconditions: ensure the ClassInfo has been cloned.
	BOOST_CHECK(a0x1->getClassInfo() != a_x->getClassInfo());
	BOOST_CHECK_EQUAL(*a0x1->getClassInfo(), *a_x->getClassInfo());
}

/**
 * Check ClassInfo attachment with append.
 */
BOOST_AUTO_TEST_CASE(testScalarAppend)
{
	Scalar *a0y = new Scalar("y");
	a0->append(a0y);

	BOOST_CHECK_EQUAL(a0y->getClassInfo(), a_y->getClassInfo());
}

BOOST_AUTO_TEST_SUITE_END()

}
