// =============================================================================
// ValidationVisitorTest.cpp
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

#include "CCfg/Attribute.h"
#include "CCfg/Object.h"
#include "CCfg/Scalar.h"
#include "CCfg/ValidationVisitor.h"
#include "CCfg/Interval.h"
#include <boost/test/auto_unit_test.hpp>
#include <sstream>

namespace CCfg
{

struct TestValidationVisitor : public AbstractValidationVisitor
{
	// typedefs
	typedef std::vector<const Attribute*> BadAttributesList;

	TestValidationVisitor(int checkFlags = CheckNoDoublets | CheckRanges | CheckTypes)
		: AbstractValidationVisitor(checkFlags),
		doubletAttributes(), outofRangeAttributes(), invalidTypeAttributes()
	{
	}
	virtual ~TestValidationVisitor() {}

	void onDoublet(const Attribute& attribute)
	{
		doubletAttributes.push_back(&attribute);
	}
	void onOutOfRange(const Scalar& scalar)
	{
		outofRangeAttributes.push_back(&scalar);
	}
	void onTypeMismatch(const Attribute& attr, const Attribute& /*other*/)
	{
		invalidTypeAttributes.push_back(&attr);
	}

	BadAttributesList doubletAttributes;
	BadAttributesList outofRangeAttributes;
	BadAttributesList invalidTypeAttributes;

private:
	TestValidationVisitor(const TestValidationVisitor& rvalue);
	TestValidationVisitor& operator=(const TestValidationVisitor& rvalue);
};


///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_validator_data
{
	Object config;
	Object* a_;
	Object* a0;
	Scalar* a_x;
	Scalar* a0x;
	Scalar* a_y;
	Scalar* a0y;

	config_validator_data()
		: config(), a_(new Object("A", "*")), a0(new Object("A", "0")),
		a_x(new Scalar("x")), a0x(new Scalar("x")), a_y(new Scalar("y")),
		a0y(new Scalar("y"))
	{
		config.append(a_);
		config.append(a0);
		a_->append(a_x);
		a_->append(a_y);
		a0->append(a0x);
		a0->append(a0y);

		a_x->setRange(Interval(0,100));
		a_y->setRange(Interval(0.0, 1.0));
		*a_x = 10;
		*a0x = 20;
		*a_y = 0.1;
		*a0y = 0.75;
	}
	~config_validator_data()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_ValidationVisitor, config_validator_data)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check if the validation correct configuration works.
 */
BOOST_AUTO_TEST_CASE(testValidationVisitor1)
{
	TestValidationVisitor v;
	config.accept(v);

	BOOST_CHECK(v.doubletAttributes.empty());
	BOOST_CHECK(v.outofRangeAttributes.empty());
	BOOST_CHECK(v.invalidTypeAttributes.empty());
}

/**
 * Check if the validation of a configuration having an attribute of the wrong type works.
 */
BOOST_AUTO_TEST_CASE(testValidationVisitor2)
{
	*a_y = 0x5u;

	TestValidationVisitor v;
	config.accept(v);

	BOOST_CHECK(v.doubletAttributes.empty());
	BOOST_CHECK_EQUAL(v.outofRangeAttributes.size(), 1u);
	BOOST_CHECK_EQUAL(v.invalidTypeAttributes.size(), 1u);
	BOOST_CHECK_EQUAL(v.outofRangeAttributes[0], a_y);
	BOOST_CHECK_EQUAL(v.invalidTypeAttributes[0], a_y);
}

/**
 * Check if the validation of a configuration having an attribute out of range works.
 */
BOOST_AUTO_TEST_CASE(testValidationVisitor3)
{
	*a0y = 2.5;

	TestValidationVisitor v;
	config.accept(v);

	BOOST_CHECK(v.doubletAttributes.empty());
	BOOST_CHECK_EQUAL(v.outofRangeAttributes.size(), 1u);
	BOOST_CHECK(v.invalidTypeAttributes.empty());
	BOOST_CHECK_EQUAL(v.outofRangeAttributes[0], a0y);
}

/**
 * Check if the validation of a configuration having a doublet works.
 */
BOOST_AUTO_TEST_CASE(testValidationVisitor4)
{
	Scalar* a0x2 = new Scalar("x");
	a0->append(a0x2);
	*a0x2 = 50;

	Scalar* a0x3 = new Scalar("x");
	a0->append(a0x3);
	*a0x3 = 75;

	TestValidationVisitor v;
	config.accept(v);

	BOOST_CHECK_EQUAL(v.doubletAttributes.size(), 3u);
	BOOST_CHECK(v.outofRangeAttributes.empty());
	BOOST_CHECK(v.invalidTypeAttributes.empty());
	BOOST_CHECK_EQUAL(v.doubletAttributes[0], a0x);
	BOOST_CHECK_EQUAL(v.doubletAttributes[1], a0x2);
	BOOST_CHECK_EQUAL(v.doubletAttributes[2], a0x3);
}

BOOST_AUTO_TEST_SUITE_END()

}
