// =============================================================================
// InheritedIterator.h
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

#ifndef CCfg_View_InheritedIterator_h_INCLUDED
#define CCfg_View_InheritedIterator_h_INCLUDED

#include "InheritanceIterator.h"

namespace CCfg
{
// forward declaration
class Attribute;
class Path;

namespace View
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * Enumerates every objects that are default objects for the specified object.
 * \note
 * Default objects appears in decreasing priority order.
 *
 * author fchateau
 */
class InheritedIterator : public InheritanceIterator
{
public:
	// construction/destruction
	InheritedIterator();
	explicit InheritedIterator(Attribute& baseObject);
	InheritedIterator(Attribute& baseObject, const Path& path);
	InheritedIterator(const InheritedIterator& value);
	virtual ~InheritedIterator();
	InheritedIterator& operator =(const InheritedIterator& value);

	// iteration
	bool hasNext();

	// other methods: see base class
};

}

}

#endif // CCfg_View_InheritedIterator_h_INCLUDED
