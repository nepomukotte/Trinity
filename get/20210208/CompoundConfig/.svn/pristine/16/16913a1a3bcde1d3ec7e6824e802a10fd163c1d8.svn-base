// =============================================================================
// Path.h
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

#ifndef CCfg_Path_h_INCLUDED
#define CCfg_Path_h_INCLUDED

#include "PathNode.h"

#include <deque>

namespace CCfg
{
class Attribute;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class describes a path to a configuration attribute.
 * \note
 * Path nodes are not modifiables, iterators are in fact const_iterators.
 * author fchateau
 */
class Path
{
	std::deque<PathNode> path;

public:
	static const int LAST = -1;

	typedef std::deque<PathNode>::iterator iterator;
	typedef std::deque<PathNode>::const_iterator const_iterator;
	typedef std::deque<PathNode>::reverse_iterator reverse_iterator;
	typedef std::deque<PathNode>::const_reverse_iterator const_reverse_iterator;

	// inner classes
	template<bool> class Inserter;
	typedef Inserter<true> BackInserter;
	typedef Inserter<false> FrontInserter;

	// construction/destruction
	Path();
	Path(const std::string& path);
	Path(const_iterator start, const_iterator end);
	Path(const Path& value);
	Path& operator =(const Path& value);
	~Path();

	// comparison
	bool operator==(const Path& right) const;
	bool operator!=(const Path& right) const;

	// iteration
	iterator begin();
	iterator end();
	const_iterator begin() const;
	const_iterator end() const;
	reverse_iterator rbegin();
	reverse_iterator rend();
	const_reverse_iterator rbegin() const;
	const_reverse_iterator rend() const;

	// getters/setters
	std::string toString() const;
	void fromString(const std::string& pathString);
	const PathNode& getTail() const;
	const PathNode& getHead() const;
	PathNode& getTail();
	PathNode& getHead();
	bool empty() const;
	int size() const;

	// path construction
	BackInserter append(const std::string& name, const std::string& index = "");
	Path& append(const PathNode& newNode);
	Path& append(const Attribute& newAttr);
	Path& append(const Path& trailingPath);

	FrontInserter prepend(const std::string& name, const std::string& index = "");
	Path& prepend(const PathNode& newNode);
	Path& prepend(const Attribute& newAttr);
	Path& prepend(const Path& leadingPath);

	Path& operator<<(const PathNode& pathNode);
	Path& operator<<(const Attribute& attr);

	// path destruction
	void removeTail();
	void clear();

	// path extraction
	Path subPath(int start, int last = LAST) const;
	Path subPath(const_iterator start, const_iterator last) const;

	// other
	Path createGenericPath() const;
};

std::ostream& operator<<(std::ostream& stream, const Path& path);


/**
 * Transient object allowing to specify a list of node to insert after for a
 * call to append or prepend.
 * It enables syntax like this: path.append("A", "0")("B", "0")("x");
 */
template<bool atEnd>
class Path::Inserter
{
	Path& path;
public:
	Inserter(Path& path);
	Inserter<atEnd>& operator()(const std::string& name, const std::string& index = "");
};




///////////////////////////////////////////////////////////////////////////////
// CCfg::Path inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Unequality operator: indicates whether two paths are different or not.
 * @param right the path to compare to this instance.
 * @return true if the paths are different, false otherwise.
 */
inline bool Path::operator!=(const Path& right) const
{
	return not operator==(right);
}

/**
 * Provides an iterator on the first node of the path.
 * @return Returns an iterator on the first node of the path.
 */
inline Path::iterator Path::begin()
{
	return path.begin();
}

/**
 * Provides an iterator beyond the last node of the path.
 * @return Returns an iterator beyond the last node of the path.
 */
inline Path::iterator Path::end()
{
	return path.end();
}

/**
 * Provides an iterator on the first node of the path.
 * @return Returns an iterator on the first node of the path.
 */
inline Path::const_iterator Path::begin() const
{
	return path.begin();
}

/**
 * Provides an iterator beyond the last node of the path.
 * @return Returns an iterator beyond the last node of the path.
 */
inline Path::const_iterator Path::end() const
{
	return path.end();
}

/**
 * Provides an iterator on the last node of the path.
 * @return Returns an iterator on the last node of the path.
 */
inline Path::reverse_iterator Path::rbegin()
{
	return path.rbegin();
}

/**
 * Provides an iterator beyond the first node of the path.
 * @return Returns an iterator beyond the first node of the path.
 */
inline Path::reverse_iterator Path::rend()
{
	return path.rend();
}

/**
 * Provides an iterator on the last node of the path.
 * @return Returns an iterator on the last node of the path.
 */
inline Path::const_reverse_iterator Path::rbegin() const
{
	return path.rbegin();
}

/**
 * Provides an iterator beyond the first node of the path.
 * @return Returns an iterator beyond the first node of the path.
 */
inline Path::const_reverse_iterator Path::rend() const
{
	return path.rend();
}

/**
 * Indicates whether this path is empty or not.
 * @return Returns true if the path is empty, false otherwise.
 */
inline bool Path::empty() const
{
	return path.empty();
}

/**
 * Returns the number of nodes in the path.
 * @return Returns the number of nodes in the path.
 */
inline int Path::size() const
{
	return path.size();
}

/**
 * Returns the last node of the path.
 * The result is unpredictable if the path is empty.
 * @return the last node of the path.
 */
inline const PathNode& Path::getTail() const
{
	return path.back();
}

/**
 * Returns the last node of the path.
 * The result is unpredictable if the path is empty.
 * @return the last node of the path.
 */
inline PathNode& Path::getTail()
{
	return path.back();
}

/**
 * Returns the first node of the path.
 * The result is unpredictable if the path is empty.
 * @return the first node of the path.
 */
inline const PathNode& Path::getHead() const
{
	return path.front();
}

/**
 * Returns the first node of the path.
 * The result is unpredictable if the path is empty.
 * @return the first node of the path.
 */
inline PathNode& Path::getHead()
{
	return path.front();
}

/**
 * Add the specified node at the end of the path.
 * @param newNode the node to append.
 * @return Returns the modified path.
 */
inline Path& Path::append(const PathNode& newNode)
{
	path.push_back(newNode);
	return *this;
}


/**
 * Creates a node from the specified name and index and adds it at the end of
 * the path.
 * @param name the name of the node to append.
 * @param index the index of the node to append.
 * @return Returns an inserter object that allow to specify more nodes to
 * insert, using parentheses operator.
 * @see template<bool> class Path::Inserter
 */
inline Path::BackInserter
Path::append(const std::string& name, const std::string& index)
{
	path.push_back(PathNode(name, index));
	return Inserter<true>(*this);
}

/**
 * Add the specified node at the beginning of the path.
 * @param newNode the node to prepend.
 * @return Returns the modified path.
 */
inline Path& Path::prepend(const PathNode& newNode)
{
	path.push_front(newNode);
	return *this;
}

/**
 * Creates a node from the specified name and index and adds it at the begining
 * of the path.
 * @param name the name of the node to prepend.
 * @param index the index of the node to prepend.
 * @return Returns an inserter object that allow to specify more nodes to
 * insert, using parentheses operator. Note that each additionnal node will be
 * prepended, so they will be placed in the reverse order of their declaration.
 * Example: path.prepend("A")("B") will produce a path like B.A.path
 * @see template<bool> class Path::Inserter
 */
inline Path::FrontInserter
Path::prepend(const std::string& name, const std::string& index)
{
	path.push_front(PathNode(name, index));
	return Inserter<false>(*this);
}

/**
 * Add the specified node at the end of the path.
 * @param pathNode the node to append.
 * @return Returns the modified path.
 */
inline Path& Path::operator<<(const PathNode& pathNode)
{
	return append(pathNode);
}

/**
 * Creates a node from the specified Attribute, and adds it at the end of
 * the path.
 * @param attr the attribute to append to the path.
 * @return Returns the modified path.
 */
inline Path& Path::operator<<(const Attribute& attr)
{
	return append(attr);
}

/**
 * Removes the last node from the path.
 */
inline void Path::removeTail()
{
	path.pop_back();
}

/**
 * Clear the path of all its node.
 */
inline void Path::clear()
{
	path.clear();
}

/**
 * Creates a subpath.
 * @param start an iterator on the first path node.
 * @param last an iterator on the last path node.
 * @return Returns the newly created subpath.
 */
inline Path Path::subPath(const_iterator start, const_iterator last) const
{
	return Path(start, last);
}

///////////////////////////////////////////////////////////////////////////////
// CCfg::Path::Inserter inline methods body
///////////////////////////////////////////////////////////////////////////////

template<bool atEnd>
Path::Inserter<atEnd>::Inserter(Path& path)
	: path(path)
{
}

template<>
inline Path::Inserter<true>&
Path::Inserter<true>::operator()(const std::string& name, const std::string& index)
{
	path.append(name, index);
	return *this;
}

template<>
inline Path::Inserter<false>&
Path::Inserter<false>::operator()(const std::string& name, const std::string& index)
{
	path.prepend(name, index);
	return *this;
}

}


#endif // CCfg_Path_h_INCLUDED
