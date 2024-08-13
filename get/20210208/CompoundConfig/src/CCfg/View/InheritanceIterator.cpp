// =============================================================================
// InheritanceIterator.cpp
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

#include "InheritanceIterator.h"
#include "Search.h"
#include "../Attribute.h"

namespace CCfg
{
namespace View
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::View::InheritanceIterator class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor: creates an uninitialized inheritance view.
 */
InheritanceIterator::InheritanceIterator()
	: baseAttribute(0), currentAttribute(), nextAttribute(), path()
{
}

/**
 * Implicit constructor: creates an inheritance view on the specified
 * object, using its path.
 * @param baseAttribute the object whose inheritance will be enumerated.
 */
InheritanceIterator::InheritanceIterator(Attribute& baseAttribute)
	: baseAttribute(&baseAttribute), currentAttribute(), nextAttribute(),
	path(baseAttribute.getPath())
{
}

/**
 * Constructor: creates an inheritance view for the specified path,
 * starting at the base object.
 * @param baseAttribute the object whose inheritance will be enumerated.
 * @param path the path pattern that resulting objects must match.
 */
InheritanceIterator::InheritanceIterator(Attribute& baseAttribute, const Path& path)
	: baseAttribute(&baseAttribute), currentAttribute(), nextAttribute(),
	path(path)
{
}

/**
 * Copy constructor: duplicates an inheritance view.
 * @param value the view to copy.
 */
InheritanceIterator::InheritanceIterator(const InheritanceIterator& value)
	: baseAttribute(value.baseAttribute), currentAttribute(value.currentAttribute),
	nextAttribute(value.nextAttribute), path(value.path)
{
}

/**
 * Destructor.
 */
InheritanceIterator::~InheritanceIterator()
{
}

/**
 * Assignment operator.
 * @param value the view to copy.
 */
InheritanceIterator& InheritanceIterator::operator =(const InheritanceIterator& value)
{
	baseAttribute = value.baseAttribute;
	currentAttribute = value.currentAttribute;
	nextAttribute = value.nextAttribute;
	path = value.path;
	return *this;
}

/**
 * Indicates whether the iterator still have objects to enumerate.
 * @return Returns true if a call to next() will succeed, and false otherwise.
 */
bool InheritanceIterator::hasNext()
{
	return false;
}

/**
 * Moves to the next object.
 * The cursor move to the first object having a lesser priority than the current.
 * @return Returns a reference to the view (it allows to chain calls).
 */
InheritanceIterator& InheritanceIterator::next()
{
	if(hasNext())
	{
		currentAttribute = nextAttribute;
	}
	else
	{
		assert(false);
	}
	return *this;
}

/**
 * Moves to the next object.
 * The cursor move to the first object having a lesser priority than the current.
 * @return Returns a reference to the view (it allows to chain calls).
 */
InheritanceIterator& InheritanceIterator::operator ++()
{
	return next();
}


/**
 * Retrieves the current object in the view.
 * \note
 * The function "next()" must have been called at least once before.
 * @return Returns the current default object in the view.
 */
Attribute& InheritanceIterator::operator *()
{
	assert(currentAttribute.getAttribute() != 0);
	return *currentAttribute.getAttribute();
}

/**
 * Allows to directly call methods on the current object in the view.
 * @return Returns a pointer to the current object in the view.
 */
Attribute* InheritanceIterator::operator ->()
{
	return currentAttribute;
}

/**
 * Retrieves the current object in the view (const version).
 * @return Returns the current object in the view.
 */
const Attribute& InheritanceIterator::operator *() const
{
	assert(currentAttribute.getAttribute() != 0);
	return *currentAttribute.getAttribute();
}

/**
 * Allows to directly call methods on the current object in the view
 * (const version).
 * @return Returns a pointer to the current object in the view.
 */
const Attribute* InheritanceIterator::operator ->() const
{
	return currentAttribute.getAttribute();
}

/**
 * Retrieves the current object in the view.
 * \note
 * The function "next()" must have been called at least once before.
 * @return Returns the current object in the view.
 */
Result InheritanceIterator::get() const
{
	assert(currentAttribute.getAttribute() != 0);
	return currentAttribute;
}

}

}
