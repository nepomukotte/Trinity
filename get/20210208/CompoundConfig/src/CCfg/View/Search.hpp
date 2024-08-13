// =============================================================================
// Search.hpp
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

#ifndef CCfg_View_SEARCH_HPP_INCLUDED
#define CCfg_View_SEARCH_HPP_INCLUDED

#include "Search.h"
#include "../Object.h"
#include <string>
#include <algorithm>
#include <cassert>

using std::string;

#undef max

namespace CCfg
{
namespace View
{

////////////////////////////////////////////////////////////////////////////////
// CCfg::View::Search inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Getter of the property "maxBacktrackLevel".
 * @return Returns the maximum number of nested backtrack during search
 * execution.
 */
template<class Node>
int Search<Node>::getMaxBacktrackLevel() const
{
	return maxBacktrackLevel;
}

/**
 * Default Constructor: create a search.
 */
template<class Node>
Search<Node>::Search(FilterFunc filter)
	: filter(filter), path(), iInstance(),
	currentBacktrackLevel(0), maxBacktrackLevel(0)
{
}

/**
 * Destructor.
 */
template<class Node>
Search<Node>::~Search()
{
}

/**
 * Searches the attribute designated by the path property of the instance.
 *
 * @param fullSearchPath the absolute path (from the configuration root) to
 * the attribute to search.
 * @param start the attribute where the search begins.
 * @param startIndexInPath the index of the path node, corresponding to the
 * startup attribute.
 * @return If the search succeeds, it returns a pointer on the attribute found.
 * Otherwise, it returns 0.
 */
template<class Node>
Node* Search<Node>::search(const Path& fullSearchPath, Node& start, int startIndexInPath)
{
	path = fullSearchPath;
	currentBacktrackLevel = 0;
	maxBacktrackLevel = 0;
	// check index validity.
	assert(startIndexInPath >= 0 and startIndexInPath < path.size());

	// we assume the startup attribute is an instance.
	iInstance = path.begin() + startIndexInPath;

	// we will begin the search the node after the startup attribute.
	Path::const_iterator iPath = iInstance + 1;
	if(iPath == path.end())
	{
		// if the startup attribute is at the end of the path, there is nothing
		// more to search.
		return &start;
	}
	else
	{
		// we begin the search, allowing backtrack until the configuration root.
		return searchOrBacktrack(iPath, start, 0, false);
	}
}

template<class Node>
Node* Search<Node>::match(const Path& fullMatchPath, Node& start, int startIndexInPath)
{
	path = fullMatchPath;
	currentBacktrackLevel = 0;
	maxBacktrackLevel = 0;

	// check index validity.
	assert(startIndexInPath >= 0 and startIndexInPath < path.size());

	Path::const_iterator iPath = path.begin() + (startIndexInPath + 1);
	if(iPath == path.end())
	{
		return &start;
	}
	else
	{
		return matchOrBacktrack(iPath, start, start.begin(), 0);
	}
}

/**
 * Searches the next attribute corresponding to a path.
 * @param fullSearchPath the path pattern of the siblings to search (Indexes
 * can be substituted by the default index at each level).
 * @param currentSibling the attribute for which we want to search siblings.
 * Although not checked, its path must conform to the path pattern, otherwise
 * it will not work correctly.
 * @return If the search succeeds, it returns a pointer on the sibling found.
 * Otherwise, it returns 0.
 */
template<class Node>
Node* Search<Node>::searchNext(const Path& fullSearchPath, Node& currentSibling)
{
	path = fullSearchPath;
	currentBacktrackLevel = 0;
	maxBacktrackLevel = 0;
	// we consider that the current sibling is at the top of the path and that
	// it is an instance.
	iInstance = path.begin() + (path.size() - 1);
	// we begin the search by backtracking.
	return backtrack(path.end(), currentSibling, 0);
}

/**
 * Searches the next attribute matching the specified path pattern.
 * The match algorithm browse the configuration in the tree order, starting
 * at the specified attribute.
 * This method can be used to determine which objects are inheriting from a default
 * object.
 * @param fullMatchPath the path to match. Star (*) indexes are considered as
 * wildcards and other indexes are matched exactly.
 * @param start the attribute where the search starts. It must be an
 * attribute that already matches the path but it will not be considered as a
 * potential result.
 * @return Returns the nearest attribute after the starting attribute
 * matching the path. The result can be used as the start attribute for a
 * subsequent call.
 */
template<class Node>
Node* Search<Node>::matchNext(const Path& fullMatchPath, Node& start)
{
	path = fullMatchPath;
	currentBacktrackLevel = 0;
	maxBacktrackLevel = 0;
	// we begin the match by backtracking, because the starting attribute is at
	// the top of the path.
	return backtrackMatch(path.end(), start, 0);
}

/**
 * Recursively searches the attribute designated by the path property, starting
 * from the specified node and attribute.
 * This algorithm cannot backtrack before the startup attribute.
 * @param iPath an iterator on the next path node to search.
 * @param current the attribute from where the search begins.
 * @param btStart indicates whether we are beginning a backtracking branch.
 * @param backtrackLimit the attribute beyond which we cannot backtrack.
 * A value of 0 indicates the only backtrack limit is the root of the configuration.
 * @return If the search succeeds, it returns a pointer on the attribute found.
 * Otherwise, it returns 0.
 */
template<class Node>
Node* Search<Node>::search(Path::const_iterator iPath,
	Node& current,
	Node* backtrackLimit,
	bool btStart)
{
	// When beginning a backtrack, we must search the default object of the
	// next path node, and not the next node itself, because it is where we
	// come from.
	string index = (btStart) ? "*" : iPath->index;
	// Search the next attribute.
	NodeIterator child = current.search(iPath->name, index);
	if(child == current.end())	// not found
	{
		// if we already found instances beyond the current point
		// and we have not already tried the default object.
		if(iInstance >= iPath and index != "*")
		{
			// we can look for the default object
			index = "*";
			// Search the next default attribute.
			child = current.search(iPath->name, index);
			if(child == current.end())	// not found
			{
				return 0;
			} // else: the algorithm goes on.
		}
		else
		{
			return 0;
		}
	}

	if((*filter)(*child))
	{
		return 0;
	}
	// if the attribute we have found comes after the last instance, we update it.
	if(iInstance < iPath)
	{
		// note: by recursion ++iInstance == iPath,
		// so the following assignation is equivalent to "advance by 1 step"
		iInstance = iPath;
	}
	// if we just found the last attribute designated by the path.
	if(++iPath == path.end())
	{
		return &*child; // algorithm successful termination
	}
	else // there are more attributes to find.
	{
		return searchOrBacktrack(iPath, *child, backtrackLimit, false);
	}
}

/**
 * Recursively backtracks and search until a backtrack branch finds the
 * attribute designated by the path.
 * @param current the attribute from which we will backtrack.
 * @param backtrackLimit the attribute beyond which we cannot backtrack.
 * A value of 0 indicates the only backtrack limit is the root of the configuration.
 * @param iPath an iterator on the next path node to search.
 * @return If the search succeeds, it returns a pointer on the attribute found.
 * Otherwise, it returns 0.
 */
template<class Node>
Node* Search<Node>::backtrack(Path::const_iterator iPath, Node& current,
	Node* backtrackLimit)
{
	Node* result;
	++currentBacktrackLevel;
	maxBacktrackLevel = std::max(maxBacktrackLevel, currentBacktrackLevel);

	Node* parent = current.getParent();
	if(parent != backtrackLimit) // backtrack can only occur above the backtrackLimit.
	{
		--iPath; // move back the path iterator.
		if(current.getIndex() == "*")
		{
			// if the index of the current attribute is *, we have depleted
			// all possibilities at this level, we must backtrack.
			result = backtrack(iPath, *parent, backtrackLimit);
		} else
		{
			// otherwise, we can start a backtrack branch with the * index.
			// the value "true" for parameter btStart will indicates to search
			// the default object.
			result = searchOrBacktrack(iPath, *parent, backtrackLimit, true);
		}
	} else
	{
		result = 0;
	}

	--currentBacktrackLevel;
	return result;
}

/**
 * Recursively searches the attribute designated by the path property, starting
 * from the specified node and attribute.
 * This algorithm can backtrack before the startup attribute if it cannot find
 * the next attribute.
 * @param iPath an iterator on the next path node to search.
 * @param current the attribute from where the search begins.
 * @param btStart indicates whether we are beginning a backtracking branch.
 * @param backtrackLimit the attribute beyond which we cannot backtrack.
 * A value of 0 indicates the only backtrack limit is the root of the configuration.
 * @return If the search succeeds, it returns a pointer on the attribute found.
 * Otherwise, it returns 0.
 */
template<class Node>
Node* Search<Node>::searchOrBacktrack(const Path::const_iterator& iPath,
	Node& current,
	Node* backtrackLimit,
	bool btStart)
{
	// nested backtracks cannot go past the current object
	Node* result = search(iPath, current, current.getParent(), btStart);
	if(result == 0)
	{
		result = backtrack(iPath, current, backtrackLimit);
	}
	return result;
}



/**
 * Recursively searches the next attribute matching the path property.
 * @param iPath an iterator on the next path node to match.
 * @param current the attribute from where the search begins.
 * @param nextChild an iterator on a child of "current". The method will start
 * its match evaluation from this attribute. It doesn't have to be a potential
 * result, any child is suitable.
 * @param backtrackLimit the attribute beyond which we cannot backtrack.
 * A value of 0 indicates the only backtrack limit is the root of the configuration.
 * @return If the search succeeds, it returns a pointer on the attribute found.
 * Otherwise, it returns 0.
 */
template<class Node>
Node* Search<Node>::match(Path::const_iterator iPath,
	Node& current,
	NodeIterator nextChild,
	Node* backtrackLimit)
{
	// enumerate every child of the current attribute
	NodeIterator iChild = nextChild;
	NodeIterator endChild = current.end();
	for(; iChild != endChild; ++iChild)
	{
		Node& child = *iChild;
		// check if the attribute name and index match those of the path node.
		if(child.getName() == iPath->name and
			(iPath->index == "*" or child.getIndex() == iPath->index))
		{
			if(not (*filter)(child))
			{
				// if we just found the last attribute designated by the path.
				if(++iPath == path.end())
				{
					return &child; // algorithm successful termination
				}
				else // there are more attributes to find.
				{
					return matchOrBacktrack(iPath, child, child.begin(), backtrackLimit);
				}
			}
		}
	}
	// no child matches the path node.
	return 0;
}

/**
 * Recursively backtracks and search until a backtrack branch finds an
 * attribute matching the path property.
 * @param current the attribute from which we will backtrack.
 * @param backtrackLimit the attribute beyond which we cannot backtrack.
 * A value of 0 indicates the only backtrack limit is the root of the configuration.
 * @param iPath an iterator on the next path node to search.
 * @return If the search succeeds, it returns a pointer on the attribute found.
 */
template<class Node>
Node* Search<Node>::backtrackMatch(Path::const_iterator iPath,
	Node& current,
	Node* backtrackLimit)
{
	Node* result;
	++currentBacktrackLevel;
	maxBacktrackLevel = std::max(maxBacktrackLevel, currentBacktrackLevel);

	Node* parent = current.getParent();
	if(parent != backtrackLimit) // backtrack can only occur above the backtrackLimit.
	{
		--iPath; // move back the path iterator.
		NodeIterator currentChild = current.getIterator();
		result = matchOrBacktrack(iPath, *parent, ++currentChild, backtrackLimit);
	}
	else
	{
		result = 0;
	}

	--currentBacktrackLevel;
	return result;
}

/**
 * Recursively searches the next attribute matching the path property.
 * This algorithm can backtrack before the startup attribute if it cannot find
 * the next attribute.
 * @param iPath an iterator on the next path node to match.
 * @param current the attribute from where the search begins.
 * @param nextChild an iterator on a child of "current". The method will start
 * its match evaluation from this attribute. It doesn't have to be a potential
 * result, any child is suitable.
 * @param backtrackLimit the attribute beyond which we cannot backtrack.
 * A value of 0 indicates the only backtrack limit is the root of the configuration.
 * @return If the search succeeds, it returns a pointer on the attribute found.
 * Otherwise, it returns 0.
 */
template<class Node>
Node* Search<Node>::matchOrBacktrack(const Path::const_iterator& iPath,
	Node& current,
	NodeIterator nextChild,
	Node* backtrackLimit)
{
	// nested backtracks cannot go past the current object
	Node* result = match(iPath, current, nextChild, current.getParent());
	if(result == 0)
	{
		return backtrackMatch(iPath, current, backtrackLimit);
	}
	else
	{
		return result;
	}
}

template<class Node>
bool Search<Node>::defaultFilter(Node& /*node*/)
{
	return false;
}


}
}

#endif // CCfg_View_SEARCH_HPP_INCLUDED
