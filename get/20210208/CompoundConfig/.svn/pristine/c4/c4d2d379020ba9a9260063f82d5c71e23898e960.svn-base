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

#include "PathNode.h"
#include "Attribute.h"
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_closure.hpp>
#include <complex>
#include <boost/spirit/include/phoenix1.hpp>
#include <sstream>

using phoenix::arg1;
using phoenix::arg2;
using phoenix::var;
using namespace boost::spirit::classic;

namespace CCfg
{

/**
 * Stream insertion operator: prints a path node on a stream.
 * @param stream the stream where the node will be print.
 * @param node the node to print.
 * @return Returns the stream in parameter, to chain calls.
 */
std::ostream& operator<<(std::ostream& stream, const PathNode& node)
{
	stream << node.name;
	if(not node.index.empty())
	{
		stream << '[' << node.index << ']';
	}
	return stream;
}

///////////////////////////////////////////////////////////////////////////////
// NodeParser class methods body
///////////////////////////////////////////////////////////////////////////////

struct PathNodeClosure : closure<PathNodeClosure, PathNode>
{
    member1 node;
};

/**
 * This class encapsulate the grammar for parsing a path node.
 */
struct NodeParser : public grammar<NodeParser, PathNodeClosure::context_t>
{
    template <typename ScannerT>
    struct definition
    {
        definition(NodeParser const& self)
        {
        	using phoenix::bind;
        	using phoenix::construct_;

        	node_parser =
			(
				node_expr =
						name
						[
						 	 phoenix::bind(&PathNode::setName)
						 	 	 (self.node, construct_<std::string>(arg1, arg2))
						]
					>> !('['
					>>	index
						[
						 	 phoenix::bind(&PathNode::setIndex)
						 	 	 (self.node, construct_<std::string>(arg1, arg2))
						]
					>>	']'),
				name = (alpha_p | '_') >> *(alpha_p | digit_p | '-' | '_'),
				index = ch_p('*') | int_p | name
			);
        }

		subrule<0>  node_expr;
		subrule<1>  name;
		subrule<2>  index;

        typedef rule<ScannerT, PathNodeClosure::context_t> rule_t;

        rule_t node_parser;
        rule_t const& start() const { return node_parser; }
    };
};

///////////////////////////////////////////////////////////////////////////////
// PathNode class methods body
///////////////////////////////////////////////////////////////////////////////


/**
 * Default constructor.
 */
PathNode::PathNode()
	: name(), index()
{
}

/**
 * Implicit constructor: initializes this node with informations from a
 * Attribute.
 * @param attribute the Attribute from which this node will bear the name.
 */
PathNode::PathNode(const Attribute& attribute)
	: name(attribute.getName()), index(attribute.getIndex())
{
}

/**
 * Copy constructor.
 * @param value the node to copy.
 */
PathNode::PathNode(const PathNode& value)
	: name(value.name), index(value.index)
{
}

/**
 * Assignment operator.
 * @param value the node to copy.
 */
PathNode& PathNode::operator =(const PathNode& value)
{
	name = value.name;
	index = value.index;
	return *this;
}

/**
 * Destructor.
 */
PathNode::~PathNode()
{
}

/**
 * Constructor: initializes an node with the specified name and index.
 * @param name the name of the attribute this node represents.
 * @param index the index of the attribute this node represents.
 */
PathNode::PathNode(const std::string& name, const std::string& index)
	: name(name), index(index)
{
}

/**
 * Constructor: initializes an node with the specified name and index.
 * @param name the name of the attribute this node represents.
 * @param index the index of the attribute this node represents.
 */
PathNode::PathNode(const std::string& name, int index)
	: name(name)
{
	std::ostringstream indexStr;
	indexStr << index;
	this->index = indexStr.str();
}

/**
 * Converts this path node to a string.
 * @return Returns a string representation of this PathNode.
 */
std::string PathNode::toString() const
{
	std::ostringstream stream;
	stream << *this;
	return stream.str();
}

/**
 * Initialize this path node by parsing a string representing a path node.
 * @param nodeStr a string of the form: name [ '[' index ']' ].
 * The name must begin by a letter or an underscore and go on with letters,
 * digits, underscores or dashes.
 * The index must be either an integer or a string having the same
 * restrictions as the name.
 * @return Returns true if the parsing was successful, false otherwise.
 * In case of unsuccessful parsing, the path node has most likely been altered.
 */
bool PathNode::fromString(const std::string& nodeStr)
{
	index.clear();
	NodeParser node;
	return parse(nodeStr.c_str(), node[var(*this) = arg1]).full;
}

}

