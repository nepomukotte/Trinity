// =============================================================================
// ClassIterator.h
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

#ifndef CCfg_View_ClassIterator_h_INCLUDED
#define CCfg_View_ClassIterator_h_INCLUDED

#include "InheritanceIterator.h"

namespace CCfg
{
// forward declarations
class Path;
class Attribute;

namespace View
{

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 *
 * @author fchateau
 */
class ClassIterator : public InheritanceIterator
{
	int attributeDepth;
public:
	ClassIterator();
	explicit ClassIterator(Attribute& someAttribute);
	ClassIterator(Attribute& root, const Path& path);
	ClassIterator(const ClassIterator& right);
	virtual ~ClassIterator();
	ClassIterator& operator=(const ClassIterator& right);

	// iteration
	bool hasNext();

	// other methods: see base class
};

}
}

#endif // CCfg_View_ClassIterator_h_INCLUDED
