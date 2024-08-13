// =============================================================================
// PathTest.cpp
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
#include "CCfg/Path.h"
#include "CCfg/Scalar.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_path_data
{
	Path path;

	config_path_data()
	{
		path.append("cfg")("A", "0")("B", "1")("x");
	}

	~config_path_data()
	{

	}
};


BOOST_FIXTURE_TEST_SUITE(CCfg_Path, config_path_data)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////


/**
 * Checks path construction and iteration.
 */
BOOST_AUTO_TEST_CASE(testIteration)
{
	Path::const_iterator iPath = path.begin();
	Path::const_iterator pEnd = path.end();
	BOOST_CHECK(iPath != pEnd);
	BOOST_CHECK_EQUAL(iPath->name, "cfg");
	BOOST_CHECK_EQUAL(iPath->index, "");
	++iPath;
	BOOST_CHECK(iPath != pEnd);
	BOOST_CHECK_EQUAL(iPath->name, "A");
	BOOST_CHECK_EQUAL(iPath->index, "0");
	++iPath;
	BOOST_CHECK(iPath != pEnd);
	BOOST_CHECK_EQUAL(iPath->name, "B");
	BOOST_CHECK_EQUAL(iPath->index, "1");
	++iPath;
	BOOST_CHECK(iPath != pEnd);
	BOOST_CHECK_EQUAL(iPath->name, "x");
	BOOST_CHECK_EQUAL(iPath->index, "");
	++iPath;
	BOOST_CHECK(iPath == pEnd);
}

/**
 * Checks path comparison
 */
BOOST_AUTO_TEST_CASE(testComparison)
{
	Path p2;
	p2.append("cfg")("A", "0")("B", "1")("x");
	BOOST_CHECK_MESSAGE(p2 == path, "checks equal paths");

	Path p3;
	p3.append("cfg")("A", "0")("B", "1");
	BOOST_CHECK_MESSAGE(p3 != path, "checks paths of different length");

	Path p4;
	p4.append("cfg")("A", "0")("B", "10")("x");
	BOOST_CHECK_MESSAGE(p4 != path, "checks paths with different index");

	Path p5;
	p5.append("cfg")("A", "0")("Z", "1")("x");
	BOOST_CHECK_MESSAGE(p5 != path, "checks paths with different name");
}

/**
 * Checks path copy construction.
 */
BOOST_AUTO_TEST_CASE(testCopyConstruction)
{
	Path p(path);
	BOOST_CHECK(p == path);
}

/**
 * Checks path assignment.
 */
BOOST_AUTO_TEST_CASE(testPathAssignment)
{
	Path p;
	// create a dummy path
	p << PathNode("Z", 4) << PathNode("x");
	//reassign it
	p = path;
	BOOST_CHECK(p == path);
}

/**
 * Checks path reverse iterator.
 */
BOOST_AUTO_TEST_CASE(testReverseIteration)
{
	Path::const_reverse_iterator iPath = path.rbegin();
	Path::const_reverse_iterator pEnd = path.rend();
	BOOST_CHECK(iPath != pEnd);
	BOOST_CHECK_EQUAL(iPath->name, "x");
	BOOST_CHECK_EQUAL(iPath->index, "");
	++iPath;
	BOOST_CHECK(iPath != pEnd);
	BOOST_CHECK_EQUAL(iPath->name, "B");
	BOOST_CHECK_EQUAL(iPath->index, "1");
	++iPath;
	BOOST_CHECK(iPath != pEnd);
	BOOST_CHECK_EQUAL(iPath->name, "A");
	BOOST_CHECK_EQUAL(iPath->index, "0");
	++iPath;
	BOOST_CHECK(iPath != pEnd);
	BOOST_CHECK_EQUAL(iPath->name, "cfg");
	BOOST_CHECK_EQUAL(iPath->index, "");
	++iPath;
	BOOST_CHECK(iPath == pEnd);
}

/**
 * Checks string conversion (and indirectly stream insertion operators)
 */
BOOST_AUTO_TEST_CASE(testToString)
{
	BOOST_CHECK_EQUAL(path.toString(), "cfg.A[0].B[1].x");
}

/**
 * Checks various getters.
 */
BOOST_AUTO_TEST_CASE(testPathGetters)
{
	BOOST_CHECK_EQUAL(path.empty(), false);
	BOOST_CHECK_EQUAL(path.size(), 4);
	BOOST_CHECK(path.getTail() == PathNode("x"));
}

/**
 * Checks append.
 */
BOOST_AUTO_TEST_CASE(testNodeAppend1)
{
	path.append(PathNode("y"));
	BOOST_CHECK(path.getTail() == PathNode("y"));
}

/**
 * Checks append.
 */
BOOST_AUTO_TEST_CASE(testNodeAppend2)
{
	Scalar y("y", "0");
	path.append(y);
	BOOST_CHECK(path.getTail() == PathNode("y", "0"));
}


/**
 * Checks append.
 */
BOOST_AUTO_TEST_CASE(testNodeAppend3)
{
	path.append("y", "1");
	BOOST_CHECK(path.getTail() == PathNode("y", "1"));
}

/**
 * Checks append.
 */
BOOST_AUTO_TEST_CASE(testNodeAppend4)
{
	Path tail;
	tail.append("y", "1")("z", "2");
	path.append(tail);
	BOOST_CHECK_EQUAL(path.toString(), "cfg.A[0].B[1].x.y[1].z[2]");
}


/**
 * Checks removeTail.
 */
BOOST_AUTO_TEST_CASE(testRemoveTail)
{
	path.removeTail();
	BOOST_CHECK(path.getTail() == PathNode("B", "1"));
}

/**
 * Checks clear.
 */
BOOST_AUTO_TEST_CASE(testPathClear)
{
	path.clear();
	BOOST_CHECK(path.empty());
	BOOST_CHECK_EQUAL(path.size(), 0);
}

/**
 * Checks subPath.
 */
BOOST_AUTO_TEST_CASE(testSubPath)
{
	Path p(path.subPath(0));
	BOOST_CHECK(p == path);

	Path p2 = path.subPath(1, 2);
	Path p3;
	p3 << PathNode("A", 0) << PathNode("B", 1);
	BOOST_CHECK(p2 == p3);
}

/**
 * Checks prepend.
 */
BOOST_AUTO_TEST_CASE(testPrepend1)
{
	path.prepend(PathNode("y"));
	BOOST_CHECK(path.getHead() == PathNode("y"));
}

/**
 * Checks prepend.
 */
BOOST_AUTO_TEST_CASE(testPrepend2)
{
	Scalar y("y", "0");
	path.prepend(y);
	BOOST_CHECK(path.getHead() == PathNode("y", "0"));
}

/**
 * Checks prepend.
 */
BOOST_AUTO_TEST_CASE(testPrepend3)
{
	path.prepend("module", "readout")("system", "daq");
	BOOST_CHECK_EQUAL(path.toString(), "system[daq].module[readout].cfg.A[0].B[1].x");
}

/**
 * Checks prepend.
 */
BOOST_AUTO_TEST_CASE(testPrepend4)
{
	Path head;
	head.append("system", "daq")("module", "readout");
	path.prepend(head);
	BOOST_CHECK_EQUAL(path.toString(), "system[daq].module[readout].cfg.A[0].B[1].x");
}

/**
 * Checks fromString
 */
BOOST_AUTO_TEST_CASE(testPathFromString)
{
	Path p;
	p.fromString("cfg.A[0].B[1].x");
	BOOST_CHECK_MESSAGE(p == path, "normal path");

	p.fromString("");
	BOOST_CHECK_MESSAGE(p.empty(), "empty path");
}

BOOST_AUTO_TEST_SUITE_END()

}
