// =============================================================================
// InheritersIterator.cpp
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

#include "InheritersIterator.h"
#include "Search.h"
#include "../Attribute.h"

namespace CCfg
{
namespace View
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::View::InheritersIterator class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
InheritersIterator::InheritersIterator()
{
}

/**
 * Implicit constructor: creates an inheritance view on the specified
 * object, using its path.
 * @param baseObject the object whose inheritance will be enumerated.
 */
InheritersIterator::InheritersIterator(Attribute& baseObject)
	: InheritanceIterator(baseObject)
{
}

/**
 * Constructor: creates an inheritance view for the specified path,
 * starting at the base object.
 * @param baseObject the object whose inheritance will be enumerated.
 * @param path the path pattern that resulting objects must match.
 */
InheritersIterator::InheritersIterator(Attribute& baseObject, const Path& path)
	: InheritanceIterator(baseObject, path)
{
}

/**
 * Copy constructor.
 */
InheritersIterator::InheritersIterator(const InheritersIterator& value)
	: InheritanceIterator(value)
{
}

/**
 * Destructor.
 */
InheritersIterator::~InheritersIterator()
{
}

/*
 *
 */
bool InheritersIterator::hasNext()
{
	Attribute* cur =
		(currentAttribute.getAttribute() == 0) ?
			baseAttribute :
			currentAttribute.getAttribute();

	if(nextAttribute.getAttribute() == currentAttribute.getAttribute())
	{
		Search<Attribute> s;
		Attribute* next = s.matchNext(path, *cur);
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
