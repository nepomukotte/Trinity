// =============================================================================
// ClassIterator.cpp
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

#include "ClassIterator.h"
#include "Search.h"
#include "../Object.h"

namespace CCfg
{
namespace View
{

////////////////////////////////////////////////////////////////////////////////
// CCfg::View::ClassIterator class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ClassIterator::ClassIterator()
	: InheritanceIterator(), attributeDepth()
{
}

/**
 * Full constructor.
 * @param someAttribute an instance of the class we want to enumerate.
 */
ClassIterator::ClassIterator(Attribute& someAttribute)
	: InheritanceIterator(*someAttribute.getRootConfig(),
		someAttribute.getPath().createGenericPath()), attributeDepth(0)
{
}

/**
 * Full constructor.
 * @param root the root of the configuration of the class to enumerate.
 * @param path the path of an instance of the class to enumerate.
 */
ClassIterator::ClassIterator(Attribute& root, const Path& path)
	: InheritanceIterator(root, path.createGenericPath()), attributeDepth(0)
{
}

/**
 * Copy constructor.
 */
ClassIterator::ClassIterator(const ClassIterator& right)
	: InheritanceIterator(right), attributeDepth(right.attributeDepth)
{
}

/**
 * Destructor.
 */
ClassIterator::~ClassIterator()
{
}

/**
 * Assignment operator.
 */
ClassIterator& ClassIterator::operator=(const ClassIterator& right)
{
	InheritanceIterator::operator=(right);
	attributeDepth = right.attributeDepth;
	return *this;
}

/*
 *
 */
bool ClassIterator::hasNext()
{
	if(nextAttribute.getAttribute() == currentAttribute.getAttribute())
	{
		Search<Attribute> s;
		Attribute* next;
		if(currentAttribute.getAttribute() == 0)
		{
			next = s.match(path, *baseAttribute, attributeDepth);
		}
		else
		{
			next = s.matchNext(path, *currentAttribute.getAttribute());
		}
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
