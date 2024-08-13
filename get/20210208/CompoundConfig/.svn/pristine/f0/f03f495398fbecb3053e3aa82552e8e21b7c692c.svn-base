// =============================================================================
// InheritedIterator.cpp
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

#include "InheritedIterator.h"
#include "Search.h"
#include "../Attribute.h"

#include <assert.h>

namespace CCfg
{
namespace View
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::View::InheritedIterator class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor: creates an uninitialized default objects view.
 */
InheritedIterator::InheritedIterator()
	: InheritanceIterator()
{
}

/**
 * Implicit constructor: creates a default objects view on the specified
 * object, using its path.
 * @param baseObject the object whose default objects will be enumerated.
 */
InheritedIterator::InheritedIterator(Attribute& baseObject)
	: InheritanceIterator(baseObject)
{
}

/**
 * Constructor: creates a default objects view for the specified path,
 * starting at the base object.
 * The base object is generally already a default object for the attribute
 * designated by the path. The attribute designated by the path may as well
 * not exist.
 * @param baseObject the object whose default objects will be enumerated.
 * @param path the path pattern that default object must match.
 */
InheritedIterator::InheritedIterator(Attribute& baseObject, const Path& path)
	: InheritanceIterator(baseObject, path)
{
}

/**
 * Copy constructor: duplicates a default objects view.
 * @param value the view to copy.
 */
InheritedIterator::InheritedIterator(const InheritedIterator& value)
	: InheritanceIterator(value)
{
}

/**
 * Destructor.
 */
InheritedIterator::~InheritedIterator()
{
}

/**
 * Assignment operator.
 * @param value the view to copy.
 */
InheritedIterator& InheritedIterator::operator =(const InheritedIterator& value)
{
	InheritanceIterator::operator=(value);
	return *this;
}

/*
 *
 */
bool InheritedIterator::hasNext()
{
	Attribute* cur =
		(currentAttribute.getAttribute() == 0) ?
			baseAttribute :
			currentAttribute.getAttribute();

	if(nextAttribute.getAttribute() == currentAttribute.getAttribute())
	{
		Search<Attribute> s;
		Attribute *next = s.searchNext(path, *cur);
		if(next == 0)
		{
			return false;
		}
		nextAttribute = Result(next, s.getMaxBacktrackLevel());
	}

	return true;
}


}

}
