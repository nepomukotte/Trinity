// =============================================================================
// StructureNode.cpp
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

#include "StructureNode.h"
#include "Attribute.h"
#include "Range.h"

namespace CCfg
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::StructureNode class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
StructureNode::StructureNode() : name(), path(), type(OBJECT), unit(), range(),
	parent(0), children()
{
}

/**
 * Browses the structure tree and calls visitor methods for each node.
 */
bool StructureNode::accept(Visitor& visitor) const
{
	if(visitor.onClassNodeEnter(*this))
	{
		for(List::const_iterator it = children.begin(), end = children.end(); it != end; ++it)
		{
			if(not it->accept(visitor))
			{
				break;
			}
		}
	}
	return visitor.onClassNodeExit(*this);
}

/**
 * Adds a node at the end of the children list of this one.
 */
void StructureNode::append(StructureNode& node)
{
	children.push_back(node);
	node.parent = this;
}

/**
 * Returns a node with the specified name in the list of children structure nodes.
 * If none can be found, it returns null.
 */
StructureNode*
StructureNode::searchChildWithName(const std::string& name) const
{
	for(List::const_iterator it = children.begin(); it != children.end(); ++it)
	{
		const StructureNode& node = *it;
		if(node.name == name)
		{
			return const_cast<StructureNode*>(&node);
		}
	}
	return 0;
}

/**
 * Generates a structure tree for the specified configuration.
 */
std::auto_ptr<StructureNode>
StructureNode::generateStructure(const Attribute& root)
{
	std::auto_ptr<StructureNode> structRoot(new StructureNode);
	buildStructureTree(*structRoot, root, false);
	return structRoot;
}


/**
 * Recursively builds a structure tree.
 */
void
StructureNode::buildStructureTree(StructureNode& node, const Attribute& attr, bool update)
{
	if(not update)
	{
		node.name = attr.getName();
		node.path = attr.getShortClassPath().toString();
		if(attr.isObject())
		{
			node.type = OBJECT;
		}
		else
		{
			node.type = TypeInfo::getAttributeTypeId(attr.getValue().getType());
			ClassInfoPtr classInfo = attr.getClassInfo();
			if(classInfo.get() != 0)
			{
				node.unit = attr.getUnit();
				Range* range = attr.getRange();
				if(range != 0)
				{
					node.range = range->toString();
				}
			}
		}
	}
	if(attr.isObject())
	{
		for(Attribute::const_iterator it = attr.begin(); it != attr.end(); ++it)
		{
			const Attribute& childAttr = *it;
			// we look for a child with the same name, because if there is one
			// they share the same equivalence class, so we must not create
			// another structure node.
			StructureNode* otherNode = node.searchChildWithName(childAttr.getName());
			if(otherNode == 0)
			{
				node.children.push_back(StructureNode());
				StructureNode& childNode = node.children.back();
				childNode.parent = &node;
				buildStructureTree(childNode, childAttr, false);
			}
			else
			{
				buildStructureTree(*otherNode, childAttr, true);
			}
		}
	}
}

/**
 * Recursively compares 2 StructureNode trees for equality.
 * Each difference is reported throught a ComparisonMonitor object.
 * @param left the first tree.
 * @param right the second tree.
 * @param monitor the callback object that receives notifications about differences.
 */
bool
StructureNode::compare(const StructureNode& left,
		const StructureNode& right,
		ComparisonMonitor& monitor)
{
	bool isEqual = true;
	if(left.type != right.type)
	{
		monitor.onDifferentType(left, right);
		isEqual = false;
	}
	if(left.unit != right.unit)
	{
		monitor.onDifferentUnit(left, right);
		isEqual = false;
	}
	if(left.range != right.range)
	{
		monitor.onDifferentRange(left, right);
		isEqual = false;
	}
	if(not left.children.empty() or not right.children.empty())
	{
		// for each left child, search for a matching right child
		for(List::const_iterator iLeft = left.children.begin(),
			lEnd = left.children.end(); iLeft != lEnd; ++iLeft)
		{
			const StructureNode& leftChild = *iLeft;
			const StructureNode* rightChild = right.searchChildWithName(leftChild.name);
			if(rightChild == 0)
			{
				monitor.onRightClassNotFound(leftChild);
				isEqual = false;
			}
			else
			{
				isEqual = compare(leftChild, *rightChild, monitor) and isEqual;
			}
		}
		// for each right child, search for a non-matching left child to report the error.
		for(List::const_iterator iRight = right.children.begin(),
			lRight = right.children.end(); iRight != lRight; ++iRight)
		{
			const StructureNode& rightChild = *iRight;
			const StructureNode* leftChild = left.searchChildWithName(rightChild.name);
			if(leftChild == 0)
			{
				monitor.onLeftClassNotFound(rightChild);
				isEqual = false;
			}
		}
	}
	return isEqual;
}

/**
 * Recursively computes the number of nodes in a structure subtree.
 */
size_t StructureNode::computeNodeCount(const StructureNode& tree)
{
	struct CounterVisitor : Visitor
	{
		size_t count;
		CounterVisitor() : Visitor(), count(0) {}
		~CounterVisitor() {}
		bool onClassNodeEnter(const StructureNode& /*classNode*/) { ++count; return true; }
		bool onClassNodeExit(const StructureNode& /*classNode*/) { return true; }
	} v;
	tree.accept(v);
	return v.count;
}

/**
 * @fn void StructureNode::ComparisonMonitor::onLeftClassNotFound(const StructureNode& rightClass)
 * Notifies that a config class only exists in the right member of the comparison.
 * @param rightClass a reference to the existing class.
 */

/**
 * @fn void StructureNode::ComparisonMonitor::onRightClassNotFound(const StructureNode& leftClass)
 * Notifies that a config class only exists in the left member of the comparison.
 * @param leftClass a reference to the existing class.
 */

/**
 * @fn void StructureNode::ComparisonMonitor::onDifferentUnit(const StructureNode& leftClass, const StructureNode& rightClass)
 * Notifies that the unit of a class is different between the left and right members of the comparison.
 * @param leftClass the left class.
 * @param rightClass the right class.
 */

/**
 * @fn void StructureNode::ComparisonMonitor::onDifferentRange(const StructureNode& leftClass, const StructureNode& rightClass)
 * Notifies that the unit of a class is different between the left and right members of the comparison.
 * @param leftClass the left class.
 * @param rightClass the right class.
 */

/**
 * @fn void StructureNode::ComparisonMonitor::onDifferentType(const StructureNode& leftClass, const StructureNode& rightClass)
 * Notifies that the unit of a class is different between the left and right members of the comparison.
 * @param leftClass the left class.
 * @param rightClass the right class.
 */


}
