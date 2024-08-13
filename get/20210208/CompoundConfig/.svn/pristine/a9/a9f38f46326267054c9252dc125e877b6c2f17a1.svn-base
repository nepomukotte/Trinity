// =============================================================================
// StructureNode.h
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

#ifndef CCfg_StructureNode_h_INCLUDED
#define CCfg_StructureNode_h_INCLUDED

#include "TypeInfo.h"
#include <string>
#include <list>
#include <memory>

namespace CCfg
{
// forward declarations
class Attribute;

/**
 * This class describe an element in a configuration structure.
 * Each equivalence class of attribute should be described by an instance.
 *
 * @author fchateau
 */
struct StructureNode
{
	struct ComparisonMonitor;
	struct Visitor;

	typedef std::list<StructureNode> List;

	std::string name;
	std::string path;
	AttributeTypeId type;
	std::string unit;
	std::string range;
	StructureNode* parent;
	List children;

	StructureNode();
	bool accept(Visitor& visitor) const;

	void append(StructureNode& node);
	StructureNode* searchChildWithName(const std::string& name) const;

	static std::auto_ptr<StructureNode> generateStructure(const Attribute& root);
	static bool compare(const StructureNode& left, const StructureNode& right, ComparisonMonitor& monitor);
	static size_t computeNodeCount(const StructureNode& tree);

private:
	static void buildStructureTree(StructureNode& node, const Attribute& attr, bool update);
};

/**
 * This is a callback interface that must be implemented to be notified of the
 * differences between 2 configuration structures.
 */
struct StructureNode::ComparisonMonitor
{
	virtual ~ComparisonMonitor() {}
	virtual void onLeftClassNotFound(const StructureNode& rightClass) = 0;
	virtual void onRightClassNotFound(const StructureNode& leftClass) = 0;
	virtual void onDifferentUnit(const StructureNode& leftClass, const StructureNode& rightClass) = 0;
	virtual void onDifferentRange(const StructureNode& leftClass, const StructureNode& rightClass) = 0;
	virtual void onDifferentType(const StructureNode& leftClass, const StructureNode& rightClass) = 0;
};

struct StructureNode::Visitor
{
	virtual ~Visitor() {}
	virtual bool onClassNodeEnter(const StructureNode& classNode) = 0;
	virtual bool onClassNodeExit(const StructureNode& classNode) = 0;
};

}

#endif // CCfg_StructureNode_h_INCLUDED
