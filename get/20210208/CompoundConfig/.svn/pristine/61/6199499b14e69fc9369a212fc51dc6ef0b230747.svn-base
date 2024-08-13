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

#include "CCfg/Path.h"
#include "CCfg/Attribute.h"
#include "CCfg/Object.h"
#include "CCfg/Scalar.h"
#include <boost/test/auto_unit_test.hpp>
#include <cstring>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////
struct config_attribute_data
{
	std::string commentBeforeX, commentAfterX, commentY, commentBegin, commentEnd;

	Object root;
	Attribute *a, *a_foo;
	Attribute *b;
	Attribute *x;
	Attribute *y;
	Attribute *u;
	Attribute *v;
	Attribute *tab[3];

	config_attribute_data() : root()
	{
		a = new Object("A", DEFAULT_OBJ_INDEX);
		b = new Object("B");
		x = new Scalar("x");
		y = new Scalar("y");
		u = new Scalar("u");
		v = new Scalar("v");
		tab[0] = new Scalar("tab", "0");
		tab[1] = new Scalar("tab", "1");
		tab[2] = new Scalar("tab", "2");

		a_foo = new Object("A", "foo");

		root.append(a);
		root.append(a_foo);

		a->append(x);
		a->append(y);
		a->append(b);

		b->append(u);
		b->append(v);
		b->append(tab[0]);
		b->append(tab[1]);
		b->append(tab[2]);

		x->setRealValue(5.5);
		y->setStringValue("Hi");
		u->setHexValue(0xBAADF00D);
		v->setBoolValue(true);
		tab[0]->setIntValue(98);
		tab[1]->setIntValue(99);
		tab[2]->setIntValue(100);

		commentBeforeX = "this is the x attribute";
		x->setPreComment(commentBeforeX);
		x->setPostComment(commentAfterX = "nothing special");

		y->setPreComment(commentY = "yet another very useful comment");

		tab[0]->setPreComment(commentBegin = "array begin");
		tab[2]->setPostComment(commentEnd = "array end");
	}

	~config_attribute_data()
	{
		// destructor of root ensures recursive clean up
	}
};

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE(CCfg_Attribute, config_attribute_data)

/**
 * Checks whether properties are ok.
 */
BOOST_AUTO_TEST_CASE(testProperties)
{
	BOOST_CHECK_EQUAL(root.getName(), "unnamed"); // root default name
	BOOST_CHECK_EQUAL(root.getIndex(), ""); // root default index

	BOOST_CHECK_EQUAL(tab[2]->getName(), "tab"); // indexed scalar name
	BOOST_CHECK_EQUAL(tab[2]->getIndex(), "2"); // indexed scalar index

	BOOST_CHECK_EQUAL(a->getIndex(), DEFAULT_OBJ_INDEX); // test default index
	// TODO: test isDefault, etc...
}

/**
 * Checks whether parent relationships are ok.
 */
BOOST_AUTO_TEST_CASE(testParent)
{
	BOOST_CHECK_EQUAL(root.getParent(), static_cast<Object*>(0));
	BOOST_CHECK_EQUAL(a->getParent(), &root);
	BOOST_CHECK_EQUAL(b->getParent(), a);
	BOOST_CHECK_EQUAL(x->getParent(), a);
	BOOST_CHECK_EQUAL(y->getParent(), a);
	BOOST_CHECK_EQUAL(u->getParent(), b);
	BOOST_CHECK_EQUAL(v->getParent(), b);
}

/**
 * Checks whether comments works.
 */
BOOST_AUTO_TEST_CASE(testComments)
{
	BOOST_CHECK_EQUAL(x->getPreComment().getText(), commentBeforeX);
	BOOST_CHECK_EQUAL(x->getPostComment().getText(), commentAfterX);

	BOOST_CHECK_EQUAL(y->getPreComment().getText(), commentY);

	BOOST_CHECK_EQUAL(tab[0]->getPreComment().getText(), commentBegin);
	BOOST_CHECK_EQUAL(tab[2]->getPostComment().getText(), commentEnd);
}

/**
 * Checks whether attribute assignment works.
 */
BOOST_AUTO_TEST_CASE(testAttributeAssignment)
{
	Scalar val;
	Attribute& x2 = val;
	x2 = *x;
	BOOST_CHECK_EQUAL(x->getName(), x2.getName());
	BOOST_CHECK_EQUAL(x->getIndex(), x2.getIndex());
	BOOST_CHECK_EQUAL(x->getPreComment().getText(), x2.getPreComment().getText());
	BOOST_CHECK_EQUAL(x->getPostComment().getText(), x2.getPostComment().getText());
	BOOST_CHECK_EQUAL(x2.getRealValue(), 5.5);
}

/**
 * Checks whether parentheses operator and automatic cast works.
 */
BOOST_AUTO_TEST_CASE(testParentheses)
{
	double xVal = root("A", "*")("x");
	BOOST_CHECK_EQUAL(xVal, 5.5);

	const std::string& yVal = root("A", "*")("y");
	BOOST_CHECK_EQUAL(yVal, "Hi");

	unsigned uVal = root("A", "*")("B")("u");
	BOOST_CHECK_EQUAL(uVal, 0xBAADF00D);

	int tab0Val = root("A", "*")("B")("tab", 0);
	BOOST_CHECK_EQUAL(tab0Val, 98);

	int tab1Val = root("A", "*")("B")("tab", "1");
	BOOST_CHECK_EQUAL(tab1Val, 99);

	int tab2Val = root("A", "*")("B")("tab", 2);
	BOOST_CHECK_EQUAL(tab2Val, 100);
}

/**
 * Checks whether value assignment works.
 */
BOOST_AUTO_TEST_CASE(testAttributeAssignmentWithValue)
{
	root("A", "*")("x") = 2.3;
	root("A", "*")("y") = 0xDEADC0DE;
	root("A", "*")("B")("u") = 42;
	root("A", "*")("B")("v") = "Hello";

	BOOST_CHECK_EQUAL(static_cast<Real>(root("A", "*")("x")), 2.3);
	BOOST_CHECK_EQUAL(static_cast<Hex>(root("A", "*")("y")), 0xDEADC0DE);
	BOOST_CHECK_EQUAL(static_cast<Int>(root("A", "*")("B")("u")), 42);
	BOOST_CHECK_EQUAL(static_cast<const String&>(root("A", "*")("B")("v")), "Hello");
	BOOST_CHECK(not strcmp(static_cast<const char*>(root("A", "*")("B")("v")), "Hello"));
}

/**
 * Checks whether full path computation is ok.
 */
BOOST_AUTO_TEST_CASE(testFullPath)
{
	Path p = u->getPath();
	BOOST_CHECK_EQUAL(p.toString(), "unnamed.A[*].B.u");
}

/**
 * Checks whether relative path computation is ok.
 */
BOOST_AUTO_TEST_CASE(testRelativePath)
{
	Path p = tab[1]->getPath(a);
	BOOST_CHECK_EQUAL(p.toString(), "B.tab[1]");
}

/**
 * Checks whether invalid path computation is ok.
 */
BOOST_AUTO_TEST_CASE(testInvalidPath)
{
	try
	{
		Path p = tab[1]->getPath(a_foo);
		BOOST_FAIL("getPath should throw NotAnAncestor exception.");
	}
	catch(Attribute::NotAnAncestor& e)
	{
		BOOST_CHECK_EQUAL(e.getAttributePath().toString(), "unnamed.A[*].B.tab[1]");
		BOOST_CHECK_EQUAL(e.getAncestorPath().toString(), "unnamed.A[foo]");
		BOOST_CHECK_EQUAL(e.getMessage(), "Error in getPath method: the attribute unnamed.A[*].B.tab[1] is not an ancestor of unnamed.A[foo]");
	}
}

/**
 * Checks whether deep equality comparison works
 */
BOOST_AUTO_TEST_CASE(testDeepEquals)
{
	BOOST_CHECK(root.equals(root));
}

/**
 * Checks whether shallow equality comparison works
 */
BOOST_AUTO_TEST_CASE(testShallowEquals)
{
	BOOST_CHECK(a->equals(*a, IDENTICAL, SHALLOW));
}

/**
 * Checks whether cloning works.
 */
BOOST_AUTO_TEST_CASE(testClone)
{
	std::auto_ptr<Object> root2 = root.clone();
	BOOST_CHECK(root2->equals(root));
	BOOST_CHECK(root.equals(*root2));
}

/**
 * Checks whether getIncludeRoot works.
 */
BOOST_AUTO_TEST_CASE(testGetIncludeRoot)
{
	BOOST_CHECK(u->getIncludeRoot()->getIncludeUrl().empty());

	root.setIncludeUrl("foo.txt");
	BOOST_CHECK_EQUAL(u->getIncludeRoot()->getIncludeUrl(), "foo.txt");
}

/**
 * Checks whether moveTo works.
 */
BOOST_AUTO_TEST_CASE(testAttributeMoveTo)
{
	b->moveTo(root);
	BOOST_CHECK(root.search("B") == b);
	BOOST_CHECK(a->search("B") == a->end());
}

/**
 * Checks whether computeAttributeCount works.
 */
BOOST_AUTO_TEST_CASE(testComputeAttributeCount)
{
	BOOST_CHECK_EQUAL(Attribute::computeAttributeCount(root), 11u);
}

BOOST_AUTO_TEST_SUITE_END()

}
