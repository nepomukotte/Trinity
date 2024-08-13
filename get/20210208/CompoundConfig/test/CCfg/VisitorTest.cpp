// =============================================================================
// VisitorTest.cpp
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
#include "CCfg/Visitor.h"
#include <boost/test/auto_unit_test.hpp>
#include <sstream>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_visitor_data
{
	Object root;

	Object *a0;
	Object *a1;
	Scalar *x2;
	Scalar *x3;
	Scalar *x4;

	Object *a0b0;
	Object *a0b1;
	Scalar *a0y2;
	Scalar *a0y3;
	Scalar *a0y4;

	Object *a1c5;
	Object *a1c6;
	Scalar *a0z7;
	Scalar *a0z8;
	Scalar *a0z9;

	config_visitor_data() : root("root")
	{
		a0 = new Object("A", "0");
		a0b0 = new Object("B", "0");
		a0b1 = new Object("B", "1");
		a0y2 = new Scalar("y", "2");
		a0y3 = new Scalar("y", "3");
		a0y4 = new Scalar("y", "4");

		a1 = new Object("A", "1");
		a1c5 = new Object("C", "5");
		a1c6 = new Object("C", "6");
		a0z7 = new Scalar("z", "7");
		a0z8 = new Scalar("z", "8");
		a0z9 = new Scalar("z", "9");
		x2 = new Scalar("x", "2");
		x3 = new Scalar("x", "3");
		x4 = new Scalar("x", "4");

		root.append(a0);
		root.append(a1);
		root.append(x2);
		root.append(x3);
		root.append(x4);

		a0->append(a0b0);
		a0->append(a0b1);
		a0->append(a0y2);
		a0->append(a0y3);
		a0->append(a0y4);

		a1->append(a1c5);
		a1->append(a1c6);
		a1->append(a0z7);
		a1->append(a0z8);
		a1->append(a0z9);
	}

	~config_visitor_data()
	{

	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Visitor, config_visitor_data)

struct TestVisitor : public Visitor
{
	std::ostringstream sequence;

	TestVisitor()
	{
	}
	virtual ~TestVisitor()
	{
	}

	virtual bool onObjectEnter(const Object& obj)
	{
		sequence << "B" << obj.getName() << obj.getIndex() << ' ';
		return true;
	}
	virtual bool onObjectExit(const Object& obj)
	{
		sequence << "E" << obj.getName() << obj.getIndex() << ' ';
		return true;
	}
	virtual bool onScalar(const Scalar& scalar)
	{
		sequence << "S" << scalar.getName() << scalar.getIndex() << ' ';
		return true;
	}
};

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Test full visit.
 */
BOOST_AUTO_TEST_CASE(testVisit)
{
	TestVisitor visitor;
	root.accept(visitor);
	BOOST_CHECK_EQUAL(visitor.sequence.str(),
		"Broot BA0 BB0 EB0 BB1 EB1 Sy2 Sy3 Sy4 EA0 "
		"BA1 BC5 EC5 BC6 EC6 Sz7 Sz8 Sz9 EA1 "
		"Sx2 Sx3 Sx4 Eroot ");
}

/**
 * Test onObjectEnter returning false.
 */
BOOST_AUTO_TEST_CASE(testOnObjectEnter)
{
	struct TestVisitor2 : public TestVisitor
	{
		bool onObjectEnter(const Object& obj)
		{
			TestVisitor::onObjectEnter(obj);
			if(obj.getName() == "A" and obj.getIndex() == "0")	return false;
			else return true;
		}
	} visitor;
	root.accept(visitor);
	BOOST_CHECK_EQUAL(visitor.sequence.str(),
		"Broot BA0 EA0 "
		"BA1 BC5 EC5 BC6 EC6 Sz7 Sz8 Sz9 EA1 "
		"Sx2 Sx3 Sx4 Eroot ");
}

/**
 * Test onObjectExit returning false.
 */
BOOST_AUTO_TEST_CASE(testOnObjectExit)
{
	struct TestVisitor3 : public TestVisitor
	{
		bool onObjectExit(const Object& obj)
		{
			TestVisitor::onObjectExit(obj);
			if(obj.getName() == "A" and obj.getIndex() == "0")	return false;
			else return true;
		}
	} visitor;
	root.accept(visitor);
	BOOST_CHECK_EQUAL(visitor.sequence.str(),
		"Broot BA0 BB0 EB0 BB1 EB1 Sy2 Sy3 Sy4 EA0 "
		"Eroot ");
}

/**
 * Test onScalar returning false.
 */
BOOST_AUTO_TEST_CASE(testOnScalar)
{
	struct TestVisitor4 : public TestVisitor
	{
		bool onScalar(const Scalar& scalar)
		{
			TestVisitor::onScalar(scalar);
			if(scalar.getName() == "y" and scalar.getIndex() == "2")	return false;
			else return true;
		}
	} visitor;
	root.accept(visitor);
	BOOST_CHECK_EQUAL(visitor.sequence.str(),
		"Broot BA0 BB0 EB0 BB1 EB1 Sy2 EA0 "
		"BA1 BC5 EC5 BC6 EC6 Sz7 Sz8 Sz9 EA1 "
		"Sx2 Sx3 Sx4 Eroot ");
}

BOOST_AUTO_TEST_SUITE_END()

}
