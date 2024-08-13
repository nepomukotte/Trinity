// =============================================================================
// Path.cpp
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

#include "Path.h"
#include "Attribute.h"
#include <sstream>
#include <algorithm>
#include <iterator>
#include <assert.h>

using std::ostringstream;

namespace CCfg
{

/**
 * Stream insertion operator: prints a path on a stream.
 * @param stream the stream where the path will be print.
 * @param path the path to print.
 * @return Returns the stream in parameter, to chain calls.
 */
std::ostream& operator<<(std::ostream& stream, const Path& path)
{
	for(Path::const_iterator iPath = path.begin(), iEnd = path.end(); iPath != iEnd; )
	{
		stream << *iPath;
		++iPath;
		if(iPath != iEnd)
		{
			stream << '.';
		}
	}
	return stream;
}

///////////////////////////////////////////////////////////////////////////////
// Path class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor: initializes an empty path.
 */
Path::Path()
	: path()
{
}

/**
 * Implicit constructor: initializes a path by parsing a string path.
 */
Path::Path(const std::string& newPath)
	: path()
{
	fromString(newPath);
}

/**
 * Constructor: initializes a path from a path range.
 */
Path::Path(const_iterator start, const_iterator end)
	: path(start, end)
{
}

/**
 * Copy constructor.
 * @param value the Path to copy.
 */
Path::Path(const Path& value)
	: path(value.path)
{
}

/**
 * Assignment operator.
 * @param value the Path to copy.
 */
Path& Path::operator =(const Path& value)
{
	path = value.path;
	return *this;
}

/**
 * Destructor.
 */
Path::~Path()
{
}

/**
 * Equality operator: checks whether two paths are the same.
 * Two paths are identical if they have the same node count, and for each index
 * the nodes are equal.
 * @param right the path to compare to this instance.
 * @return true if the paths are the same, false otherwise.
 */
bool Path::operator==(const Path& right) const
{
	if(size() != right.size())
	{
		return false;
	}

	const_iterator iLeft = begin();
	const_iterator lEnd = end();
	const_iterator iRight = right.begin();
	const_iterator rEnd = right.end();

	for(;iLeft != lEnd and iRight != rEnd;
		++iLeft, ++iRight)
	{
		if(*iLeft != *iRight)
		{
			return false;
		}
	}
	// the test: size() == right.size() implies both paths have reached their end.
	assert(iLeft == lEnd and iRight == rEnd);
	return true;
}

/**
 * Computes the string representation of the path.
 * @return Returns the path string.
 */
std::string Path::toString() const
{
	ostringstream pathStream;
	pathStream << *this;
	return pathStream.str();
}

/**
 * Parses a string representation of a path.
 * @param pathString a string representing a path.
 */
void Path::fromString(const std::string& pathString)
{
	clear();
	std::string::const_iterator begin = pathString.begin();
	std::string::const_iterator end = begin;
	while(end != pathString.end())
	{
		end = find(begin, pathString.end(), '.');
		PathNode node;
		node.fromString(std::string(begin, end));
		append(node);
		begin = end + 1;
	}
}

/**
 * Creates a node from the specified Attribute, and adds it at the end of
 * the path.
 * @param newAttr the attribute to append to the path.
 * @return Returns the modified path.
 */
Path& Path::append(const Attribute& newAttr)
{
	path.push_back(PathNode(newAttr.getName(), newAttr.getIndex()));
	return *this;
}

/**
 * Appends a path to this one.
 * @param trailingPath the path to put as the tail of *this.
 * @return Returns the resulting path.
 */
Path& Path::append(const Path& trailingPath)
{
	std::copy(trailingPath.begin(), trailingPath.end(), std::back_inserter(path));
	return *this;
}

/**
 * Creates a node from the specified Attribute, and adds it at the
 * beginning of the path.
 * @param newAttr the attribute to prepend to the path.
 * @return Returns the modified path.
 */
Path& Path::prepend(const Attribute& newAttr)
{
	path.push_front(PathNode(newAttr.getName(), newAttr.getIndex()));
	return *this;
}

/**
 * Prepends a path to this one.
 * @param leadingPath the path to put as the head of *this.
 * @return Returns the resulting path.
 */
Path& Path::prepend(const Path& leadingPath)
{
	std::copy(leadingPath.rbegin(), leadingPath.rend(), std::front_inserter(path));
	return *this;
}

/**
 * Creates a subpath.
 * @param start the index of the first path node.
 * @param last the index of the last path node. The value Path::LAST
 * specify to use the last node.
 * @return Returns the newly created subpath.
 */
Path Path::subPath(int start, int last) const
{
	if(last == LAST)
	{
		last = path.size() - 1;
	}

	const_iterator iStart = path.begin() + start;
	const_iterator iEnd = path.begin() + last + 1;

	return Path(iStart, iEnd);
}

/**
 * Transforms a path to a generic path (that references any instances of a class).
 * A generic path has wildcard indexes.
 */
Path Path::createGenericPath() const
{
	Path genericPath = *this;
	for(Path::iterator iPath = genericPath.begin(); iPath != genericPath.end(); ++iPath)
	{
		iPath->index = "*";
	}
	return genericPath;
}


}

