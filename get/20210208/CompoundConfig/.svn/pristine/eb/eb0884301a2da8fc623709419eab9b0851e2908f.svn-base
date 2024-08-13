// =============================================================================
// StructureNodeTest.cpp
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

#include "CCfg/StructureNode.h"
#include "CCfg/Object.h"
#include "CCfg/Attribute.h"
#include "CCfg/Value.h"
#include "CCfg/Enum.h"
#include "CCfg/Interval.h"
#include <boost/test/auto_unit_test.hpp>

namespace CCfg
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct structure_node_data
{
	Object cfg1;
	structure_node_data() : cfg1("cfg1")
	{
		Attribute* attr = &cfg1.appendObject("A", "*");
		attr->appendScalar("x", "0") = 1.0;
		attr->appendScalar("x", "1")= 2.0;
		Attribute* attr2 = &attr->appendScalar("x", "2");
		attr2->setValue(Value(3.0));
		attr2->setUnit("m");
		attr2->setRange(Interval(-10, 10));
		attr2 = &attr->appendScalar("type");
		attr2->setValue(Value(2));
		Enum e;
		e = 1, 2, 3;
		attr2->setRange(e);
		attr = &cfg1.appendObject("A", "0");
		attr->appendScalar("x", "0") = 4.0;
		attr->appendScalar("x", "1") = 5.0;
		attr->appendScalar("address") = 0xABCDu;
	}
	~structure_node_data()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_StructureNode, structure_node_data)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 *
 */
BOOST_AUTO_TEST_CASE(testStructureNode1)
{
	std::auto_ptr<StructureNode> struct1 = StructureNode::generateStructure(cfg1);
	BOOST_CHECK_EQUAL(StructureNode::computeNodeCount(*struct1), 5u);
	BOOST_CHECK_EQUAL(struct1->name, "cfg1");
	BOOST_CHECK_EQUAL(struct1->path, "cfg1");
	BOOST_CHECK_EQUAL(struct1->type, OBJECT);
	BOOST_CHECK_EQUAL(struct1->unit, "");
	BOOST_CHECK_EQUAL(struct1->range, "");
	BOOST_CHECK_EQUAL(struct1->children.size(), 1u);
	StructureNode* a = &struct1->children.front();
	BOOST_CHECK_EQUAL(a->name, "A");
	BOOST_CHECK_EQUAL(a->path, "cfg1.A");
	BOOST_CHECK_EQUAL(a->type, OBJECT);
	BOOST_CHECK_EQUAL(a->unit, "");
	BOOST_CHECK_EQUAL(a->range, "");
	BOOST_CHECK_EQUAL(a->children.size(), 3u);
	StructureNode::List::iterator itAttr = a->children.begin();
	StructureNode* s = &*itAttr;
	BOOST_CHECK_EQUAL(s->name, "x");
	BOOST_CHECK_EQUAL(s->path, "cfg1.A.x");
	BOOST_CHECK_EQUAL(s->type, REAL_SCALAR);
	BOOST_CHECK_EQUAL(s->unit, "m");
	BOOST_CHECK_EQUAL(s->range, "[-10, 10]");
	BOOST_CHECK_EQUAL(s->children.size(), 0u);
	s = &*++itAttr;
	BOOST_CHECK_EQUAL(s->name, "type");
	BOOST_CHECK_EQUAL(s->path, "cfg1.A.type");
	BOOST_CHECK_EQUAL(s->type, INT_SCALAR);
	BOOST_CHECK_EQUAL(s->unit, "");
	BOOST_CHECK_EQUAL(s->range, "{1, 2, 3}");
	BOOST_CHECK_EQUAL(s->children.size(), 0u);
	s = &*++itAttr;
	BOOST_CHECK_EQUAL(s->name, "address");
	BOOST_CHECK_EQUAL(s->path, "cfg1.A.address");
	BOOST_CHECK_EQUAL(s->type, HEX_SCALAR);
	BOOST_CHECK_EQUAL(s->unit, "");
	BOOST_CHECK_EQUAL(s->range, "");
	BOOST_CHECK_EQUAL(s->children.size(), 0u);
}

/**
 *
 */
BOOST_AUTO_TEST_CASE(testStructureNode2)
{
	StructureNode root;
	root.name = "cfg1";
	root.path = "cfg1";
	root.type = OBJECT;
	StructureNode a;
	a.name = "A";
	a.path = "cfg1.A";
	a.type = OBJECT;
	root.append(a);
	StructureNode& ra = root.children.back();
	StructureNode b;
	b.name = "B";
	b.path = "cfg1.B";
	b.type = OBJECT;
	root.append(b); // difference
	StructureNode x;
	x.name = "x";
	x.path = "cfg1.A.x";
	x.type = INT_SCALAR; // difference
	x.unit = "m";
	x.range = "[-10, 10]";
	ra.append(x);
	StructureNode type;
	type.name = "type";
	type.path = "cfg1.A.type";
	type.type = INT_SCALAR;
	type.range = "{0, 1, 2, 3}"; // difference
	ra.append(type);
	// difference: no address attribute

	struct TestCompMonitor : StructureNode::ComparisonMonitor
	{
		~TestCompMonitor() {}
		void onLeftClassNotFound(const StructureNode& rightClass)
		{
			BOOST_CHECK_EQUAL(rightClass.name, "address");
		}
		void onRightClassNotFound(const StructureNode& leftClass)
		{
			BOOST_CHECK_EQUAL(leftClass.name, "B");
		}
		void onDifferentUnit(const StructureNode& /*leftClass*/, const StructureNode& /*rightClass*/)
		{
			BOOST_FAIL("units should not be different");
		}
		void onDifferentRange(const StructureNode& leftClass, const StructureNode& rightClass)
		{
			BOOST_CHECK_EQUAL(leftClass.path, rightClass.path);
			BOOST_CHECK_EQUAL(rightClass.name, "type");
			BOOST_CHECK_EQUAL(leftClass.range, "{0, 1, 2, 3}");
		}
		void onDifferentType(const StructureNode& leftClass, const StructureNode& rightClass)
		{
			BOOST_CHECK_EQUAL(leftClass.path, rightClass.path);
			BOOST_CHECK_EQUAL(leftClass.name, "x");
			BOOST_CHECK_EQUAL(leftClass.type, INT_SCALAR);
			BOOST_CHECK_EQUAL(rightClass.type, REAL_SCALAR);
		}
	};
	TestCompMonitor monitor;
	std::auto_ptr<StructureNode> root2 = StructureNode::generateStructure(cfg1);
	bool isEqual = StructureNode::compare(root, *root2, monitor);
	BOOST_CHECK(not isEqual);
	isEqual = StructureNode::compare(*root2, *root2, monitor);
	BOOST_CHECK(isEqual);
}

BOOST_AUTO_TEST_SUITE_END()

}
