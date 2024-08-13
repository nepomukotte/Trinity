// =============================================================================
// PathNode.cpp
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

#include "CCfg/PathNode.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_pathnode_data
{
	PathNode n;
	std::string nodeStr1, nodeStr2, nodeStr3, nodeStr4, nodeStr5;
	std::string nodeStr6, nodeStr7, nodeStr8;

	config_pathnode_data()
		: n(),
		nodeStr1("hello-world[0]"),
		nodeStr2("_test[foo-bar]"),
		nodeStr3("daq[*]"),
		nodeStr4("gang-of-4"),
		nodeStr5("foo[]"),
		nodeStr6("foo()"),
		nodeStr7("bar$[*]"),
		nodeStr8("bar[$42]")
	{
	}

	~config_pathnode_data()
	{
	}
};


BOOST_FIXTURE_TEST_SUITE(CCfg_PathNode, config_pathnode_data)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(testNodeFromString1)
{
	BOOST_CHECK(n.fromString(nodeStr1));
	BOOST_CHECK_EQUAL(n.name, "hello-world");
	BOOST_CHECK_EQUAL(n.index, "0");
}

BOOST_AUTO_TEST_CASE(testNodeFromString2)
{
	BOOST_CHECK(n.fromString(nodeStr2));
	BOOST_CHECK_EQUAL(n.name, "_test");
	BOOST_CHECK_EQUAL(n.index, "foo-bar");
}

BOOST_AUTO_TEST_CASE(testNodeFromString3)
{
	BOOST_CHECK(n.fromString(nodeStr3));
	BOOST_CHECK_EQUAL(n.name, "daq");
	BOOST_CHECK_EQUAL(n.index, "*");
}

BOOST_AUTO_TEST_CASE(testNodeFromString4)
{
	BOOST_CHECK(n.fromString(nodeStr4));
	BOOST_CHECK_EQUAL(n.name, "gang-of-4");
	BOOST_CHECK_EQUAL(n.index, "");
}

BOOST_AUTO_TEST_CASE(testNodeFromString5)
{
	BOOST_CHECK(not n.fromString(nodeStr5));
	BOOST_CHECK_EQUAL(n.name, "foo");
	BOOST_CHECK_EQUAL(n.index, "");
}

BOOST_AUTO_TEST_CASE(testNodeFromString6)
{
	BOOST_CHECK(not n.fromString(nodeStr6));
	BOOST_CHECK_EQUAL(n.name, "foo");
	BOOST_CHECK_EQUAL(n.index, "");
}

BOOST_AUTO_TEST_CASE(testNodeFromString7)
{
	BOOST_CHECK(not n.fromString(nodeStr7));
	BOOST_CHECK_EQUAL(n.name, "bar");
	BOOST_CHECK_EQUAL(n.index, "");
}

BOOST_AUTO_TEST_CASE(testNodeFromString8)
{
	BOOST_CHECK(not n.fromString(nodeStr8));
	BOOST_CHECK_EQUAL(n.name, "bar");
	BOOST_CHECK_EQUAL(n.index, "");
}

BOOST_AUTO_TEST_SUITE_END()

}
