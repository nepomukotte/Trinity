// =============================================================================
// Search.h
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

#ifndef CCfg_View_Search_h_INCLUDED
#define CCfg_View_Search_h_INCLUDED

#include "../Path.h"

namespace CCfg
{
namespace View
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class provide algorithms to search attributes conforming to a path,
 * using "default value" semantics.
 *
 * @author fchateau
 */
template<class Node>
class Search
{
public:
	typedef bool (*FilterFunc)(Node&);
private:
	FilterFunc filter;
	Path path;	///< the path which is the objective of the search.
	Path::const_iterator iInstance;
	/**<
	 * Transient attribute which indicates how deep is the last object instance
	 * we have encountered so far in the current search algorithm.
	 *
	 * A configuration object instance is a concrete object: not a default
	 * object. Therefore, its index is not '*'.
	 *
	 * In a search, we must find an instance for each node of the path,
	 * because a default object applies to true instances, but has no effect
	 * if there are no instances.
	 *
	 * Example: if we search A[1], which does not exist, A[*] will not replace
	 * it and should be ignored.
	 * On the other hand, if an empty A[0] exists, it will inherit A[*]
	 * attributes.
	 * A[1] could also be declared in a default object of its parent, it
	 * would be considered as a valid instance in a search for Z[0].A[1].x
	 * Z[*]
	 * 	  A[1]
	 *       x
	 * Z[0]
	 *    A[*]
	 *       x
	 *    A[0]
	 * It is worth to notice that default objects are not retroactive:
	 * Z[0].A[*] does not apply to Z[*].A[1].
	 */
	 int currentBacktrackLevel;	///< the current number of nested backtrack during algorithm execution.
	 int maxBacktrackLevel;		///< the maximum value reached by currentBacktrackLevel.

public:
	typedef typename Node::iterator NodeIterator;

	// construction/destruction
	Search(FilterFunc filter = &defaultFilter);
	~Search();

	// facade search methods
	Node* search(const Path& fullSearchPath, Node& start, int startIndexInPath);

	Node* searchNext(const Path& fullSearchPath, Node& currentSibling);

	Node* match(const Path& fullMatchPath, Node& start, int startIndexInPath);

	Node* matchNext(const Path& fullMatchPath, Node& start);

	int getMaxBacktrackLevel() const;

private:
	// non-copyable: "Search" is made to be a transient object.
	Search(const Search<Node>& value);
	Search<Node>& operator=(const Search<Node>& value);

	// search implementation methods
	Node* search(Path::const_iterator iPath, Node& current,
		Node* backtrackLimit, bool btStart);

	Node* backtrack(Path::const_iterator iPath, Node& current,
		Node* backtrackLimit);

	Node* searchOrBacktrack(const Path::const_iterator& iPath,
		Node& current, Node* backtrackLimit, bool btStart);

	// match
	Node* match(Path::const_iterator iPath, Node& current,
		NodeIterator nextChild, Node* backtrackLimit);

	Node* backtrackMatch(Path::const_iterator iPath, Node& current,
		Node* backtrackLimit);

	Node* matchOrBacktrack(const Path::const_iterator& iPath, Node& current,
		NodeIterator nextChild, Node* backtrackLimit);

	static bool defaultFilter(Node& node);
};


}

}

#endif // CCfg_View_Search_h_INCLUDED
