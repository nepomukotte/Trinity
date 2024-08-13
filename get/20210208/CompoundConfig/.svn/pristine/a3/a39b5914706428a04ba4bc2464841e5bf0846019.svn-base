// =============================================================================
// ObjectTest.cpp
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
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_object_data
{
	Object root;
	Object *a, *a_foo;
	Object *b;
	Scalar *x;
	Scalar *y;
	Scalar *u;
	Scalar *v;
	Scalar *tab[3];

	config_object_data() : root()
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
	}

	~config_object_data()
	{
		// destructor of root ensures recursive clean up
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Object, config_object_data)


///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Check whether size and empty functions works
 */
BOOST_AUTO_TEST_CASE(testSize)
{
	BOOST_CHECK_MESSAGE(root.size() == 2, "check root.size");
	BOOST_CHECK_MESSAGE(a_foo->size() == 0, "check a[foo].size");

	BOOST_CHECK_MESSAGE(root.empty() == false, "check root.empty");
	BOOST_CHECK_MESSAGE(a_foo->empty() == true, "check a[foo].empty");
}

/**
 * Checks whether forward iteration on children works.
 */
BOOST_AUTO_TEST_CASE(testForwardIteration)
{
	Attribute::iterator it = root.begin();
	Attribute::iterator end = root.end();
	BOOST_CHECK_MESSAGE(it == a, "check begin iterator");
	++it;
	BOOST_CHECK_MESSAGE(it == a_foo, "check forward iteration");
	it++;
	BOOST_CHECK_MESSAGE(it == end, "check iteration ends");
}

/**
 * Checks whether backward iteration on children works.
 */
BOOST_AUTO_TEST_CASE(testBackwardIteration)
{
	Attribute::iterator it = b->begin();
	BOOST_CHECK_MESSAGE(it == u, "check begin iterator");
	++it;
	BOOST_CHECK_MESSAGE(it == v, "check forward iteration");
	--it;
	BOOST_CHECK_MESSAGE(it == u, "check backward iteration");
}

/**
 * Checks whether searches work.
 */
BOOST_AUTO_TEST_CASE(testSearch)
{
	#define CATCH_ATTRIBUTE_NOT_FOUND(expectedMsg, expr)\
		try\
		{\
			expr;\
			BOOST_FAIL("Should throw Attribute::NotFound");\
		} catch(Attribute::NotFound& e)\
		{\
			BOOST_CHECK_EQUAL(e.getMessage(), expectedMsg); \
		}

	Attribute::iterator itA = root.find("A", DEFAULT_OBJ_INDEX);
	BOOST_CHECK_MESSAGE(itA == a, "check find A[*]");

	Attribute::iterator itA_foo = root.find("A", "foo");
	BOOST_CHECK_MESSAGE(itA_foo == a_foo, "check find A[foo]");

	Attribute::iterator itB = a->find("B");
	BOOST_CHECK_MESSAGE(itB == b, "check find B");

	CATCH_ATTRIBUTE_NOT_FOUND("Cannot find the attribute C in the object: unnamed.A[*]",
		Attribute::iterator itC = a->find("C"));

	CATCH_ATTRIBUTE_NOT_FOUND("Cannot find the attribute A[bar] in the object: unnamed",
		Attribute::iterator itA = root.find("A", "bar"));

	Attribute::iterator itTab = b->search("tab", 1);
	BOOST_CHECK_MESSAGE(itTab == tab[1], "Check search tab[1]");
	BOOST_CHECK_MESSAGE(++itTab == tab[2], "Check tab[1] iterator");

	itTab = b->search("tab", 3);
	BOOST_CHECK_MESSAGE(itTab == b->end(), "Check search tab[3] fails");
}

/**
 * Checks whether append works.
 */
BOOST_AUTO_TEST_CASE(testAttributeAppend1)
{
	Scalar *tail = new Scalar("tail");
	b->append(tail);
	Attribute::iterator iTail = b->find("tail");
	BOOST_CHECK_MESSAGE(++iTail == b->end(), "Check tail is the last attribute");
}

/**
 * Checks whether prepend works.
 */
BOOST_AUTO_TEST_CASE(testPrepend)
{
	Scalar *head = new Scalar("head");
	b->prepend(head);
	Attribute::iterator iHead = b->begin();
	BOOST_CHECK_MESSAGE(iHead == head, "Check head is the first attribute");
}

/**
 * Checks whether insert works.
 */
BOOST_AUTO_TEST_CASE(testInsert)
{
	// insert at the beginning
	Scalar *head = new Scalar("head");
	b->insert(b->begin(), head);
	Attribute::iterator iHead = b->begin();
	BOOST_CHECK_MESSAGE(iHead == head, "Check head is the first attribute");

	// insert at the end
	Scalar *tail = new Scalar("tail");
	b->insert(b->end(), tail);
	Attribute::iterator iTail = b->find("tail");
	BOOST_CHECK_MESSAGE(++iTail == b->end(), "Check tail is the last attribute");

	// insert in the middle
	Scalar *separator = new Scalar("separator");
	b->insert(b->find("tab", "0"), separator);
	Attribute::iterator iTab0 = b->find("tab", "0");
	BOOST_CHECK_MESSAGE(--iTab0 == separator, "Check separator comes before tab[0]");
}

/**
 * Checks whether remove works.
 */
BOOST_AUTO_TEST_CASE(testRemove1)
{
	a->remove(x);
	BOOST_CHECK_MESSAGE(a->search("x") == a->end(), "Check if x can no longer be found");
	a->remove("y");
	BOOST_CHECK_MESSAGE(a->search("y") == a->end(), "Check if y can no longer be found");
	Attribute::iterator it = a->begin();
	BOOST_CHECK_MESSAGE(it == b, "Check if B is now the first object");
}

/**
 * Checks whether remove works.
 */
BOOST_AUTO_TEST_CASE(testRemove2)
{
	int sizeA = a->size();
	std::auto_ptr<Attribute> b2 = a->remove(b);
	BOOST_CHECK_MESSAGE(a->search("B") == a->end(), "Check B is not a child of A anymore");
	BOOST_CHECK_EQUAL(sizeA - 1, a->size()); // Check if A.size has been decremented
	BOOST_CHECK_EQUAL(b2.get(), b); // Check if the retrieved B is the same as the original
	BOOST_CHECK_EQUAL(b2->getParent(), static_cast<Object*>(0)); // Check if B no longer has a parent
}

/**
 * Checks whether clear works.
 */
BOOST_AUTO_TEST_CASE(testObjectClear)
{
	a->clear();
	BOOST_CHECK_EQUAL(a->empty(), true); // Check if A is empty
}

/**
 * Checks whether copy construction works.
 */
BOOST_AUTO_TEST_CASE(testCopyConstruct)
{
	std::auto_ptr<Object> a2(new Object(*a));
	BOOST_CHECK_EQUAL(a->getName(), a2->getName());
	BOOST_CHECK_EQUAL(a->getIndex(), a2->getIndex());
	BOOST_CHECK_EQUAL(a->getIncludeUrl(), a2->getIncludeUrl());
	BOOST_CHECK_EQUAL(a->getPreComment().getText(), a2->getPreComment().getText());
	BOOST_CHECK_EQUAL(a->getPostComment().getText(), a2->getPostComment().getText());
	BOOST_CHECK_EQUAL(a2->empty(), true);
}

/**
 * Checks whether append works.
 */
BOOST_AUTO_TEST_CASE(testAttributeAppend2)
{
	Attribute& c = a->append(OBJECT_INST, "C");
	BOOST_CHECK(a->find("C") == &c);
	BOOST_CHECK(typeid(c) == typeid(Object));
}

/**
 * Checks whether append works.
 */
BOOST_AUTO_TEST_CASE(testAttributeAppend3)
{
	Attribute& z = a->append(SCALAR_INST, "z");
	BOOST_CHECK(a->find("z") == &z);
	BOOST_CHECK(typeid(z) == typeid(Scalar));
}

/**
 * Checks whether require works.
 */
BOOST_AUTO_TEST_CASE(testRequire1)
{
	Path p;
	p << PathNode("C") << PathNode("D");
	Attribute& d = a->require(OBJECT_INST, p);
	Attribute* c = a->find("C");
	BOOST_CHECK(c->find("D") == &d);
	BOOST_CHECK(typeid(d) == typeid(Object));
}

/**
 * Checks whether require works.
 */
BOOST_AUTO_TEST_CASE(testRequire2)
{
	Path p;
	p << PathNode("B") << PathNode("C");
	Attribute& c = a->require(OBJECT_INST, p);
	BOOST_CHECK(a->find("B") == b);
	BOOST_CHECK(b->find("C") == &c);
	BOOST_CHECK(typeid(c) == typeid(Object));
}

/**
 * Checks whether require works.
 */
BOOST_AUTO_TEST_CASE(testRequire3)
{
	Path p;
	p << PathNode("C") << PathNode("x");
	Attribute& new_x = a->require(SCALAR_INST, p);
	Attribute* c = a->find("C");
	BOOST_CHECK(c->find("x") == &new_x);
	BOOST_CHECK(typeid(new_x) == typeid(Scalar));
}

/**
 * Checks whether require works.
 */
BOOST_AUTO_TEST_CASE(testRequire4)
{
	Path p;
	p << PathNode("B") << PathNode("u");
	Attribute& new_u = a->require(SCALAR_INST, p);
	Attribute* new_b = a->find("B");
	BOOST_CHECK_EQUAL(b, new_b);
	BOOST_CHECK_EQUAL(u, &new_u);
	BOOST_CHECK(typeid(new_u) == typeid(Scalar));
}

/**
 * Check whether merge works.
 */
BOOST_AUTO_TEST_CASE(testMerge)
{
	std::auto_ptr<Attribute> attr = a->remove(b);
	a->merge(std::auto_ptr<Object>(dynamic_cast<Object*>(attr.release())));
	BOOST_CHECK(a->search("u") == u);
	BOOST_CHECK(a->search("v") == v);
	BOOST_CHECK(a->search("B") == a->end());
	BOOST_CHECK_EQUAL(a->size(), 7);
}

/**
 * Check whether import works.
 */
BOOST_AUTO_TEST_CASE(testImport)
{
	std::auto_ptr<Attribute> attr = a->remove(b);
	a->import(std::auto_ptr<Object>(dynamic_cast<Object*>(attr.release())));
	BOOST_CHECK(a->search("x") == a->end());
	BOOST_CHECK(a->search("u") == u);
	BOOST_CHECK(a->search("v") == v);
	BOOST_CHECK(a->search("B") == a->end());
	BOOST_CHECK_EQUAL(a->size(), 5);
}

BOOST_AUTO_TEST_SUITE_END()

}
