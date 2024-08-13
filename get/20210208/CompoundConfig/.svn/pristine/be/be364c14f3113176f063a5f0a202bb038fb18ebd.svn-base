// =============================================================================
// PathNode.h
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

#ifndef CCfg_PathNode_h_INCLUDED
#define CCfg_PathNode_h_INCLUDED

#include <string>

namespace CCfg
{
class Attribute;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * Describes the characteristics of a path node (corresponding to a
 * CCfg::Attribute).
 * author fchateau
 */
struct PathNode
{
	std::string name;	///< the name of the attribute represented by this node
	std::string index;	///< the index of the attribute represented by this node

	PathNode();
	explicit PathNode(const Attribute& attribute);
	explicit PathNode(const std::string& name, const std::string& index = "");
	PathNode(const std::string& name, int index);

	PathNode(const PathNode& value);
	PathNode& operator =(const PathNode& value);
	~PathNode();

	void setName(const std::string& newName) { name = newName; }
	void setIndex(const std::string& newIndex) { index = newIndex; }

	std::string toString() const;
	bool fromString(const std::string& nodeStr);

	bool operator==(const PathNode& right) const;
	bool operator!=(const PathNode& right) const;
	bool operator<(const PathNode& right) const;
};

std::ostream& operator<<(std::ostream& stream, const PathNode& node);

///////////////////////////////////////////////////////////////////////////////
// inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Equality operator: indicates whether two path nodes are the same.
 * @param right the node to compare to this instance.
 * @return Returns true if the node are equivalent, false otherwise.
 */
inline bool PathNode::operator==(const PathNode& right) const
{
	return name == right.name and index == right.index;
}

/**
 * Unequality operator: indicates whether two path nodes are different.
 * @param right the node to compare to this instance.
 * @return Returns true if the node are different, false otherwise.
 */
inline bool PathNode::operator!=(const PathNode& right) const
{
	return name != right.name or index != right.index;
}

/**
 * Ordering operator: defines the total order in a sequence of PathNode.
 * @param right the node to compare to this instance.
 * @return Returns true if the left node is lexicographically inferior
 * to the right node, false otherwise.
 */
inline bool PathNode::operator<(const PathNode& right) const
{
	return name < right.name or (name == right.name and index < right.index);
}

}


#endif // CCfg_PathNode_h_INCLUDED
