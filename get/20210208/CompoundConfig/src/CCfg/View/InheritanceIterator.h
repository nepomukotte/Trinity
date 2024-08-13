// =============================================================================
// InheritanceIterator.h
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

#ifndef CCfg_View_InheritanceIterator_h_INCLUDED
#define CCfg_View_InheritanceIterator_h_INCLUDED

#include "Result.h"
#include "../Path.h"

namespace CCfg
{
namespace View
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * A base class for classes that enumerate objects having an inheritance
 * relationship with a specified object.
 *
 * @author fchateau
 */
class InheritanceIterator
{
protected:
	Attribute* baseAttribute;	///< the attribute which is examined.
	Result currentAttribute;	///< the current attribute in iteration.
	Result nextAttribute;		///< the next attribute in iteration.
	Path path;				///< the path pattern to inheriting/inherited objects.

public:
	// construction/destruction
	InheritanceIterator();
	explicit InheritanceIterator(Attribute& baseAttribute);
	InheritanceIterator(Attribute& baseAttribute, const Path& path);
	InheritanceIterator(const InheritanceIterator& view);
	virtual ~InheritanceIterator();
	InheritanceIterator& operator=(const InheritanceIterator& view);

	// iteration
	virtual bool hasNext() = 0;
	InheritanceIterator& next();
	InheritanceIterator& operator ++();

	// dereference
	Attribute& operator *();
	Attribute* operator ->();
	const Attribute& operator *() const;
	const Attribute* operator ->() const;
	Result get() const;
};

}

}

#endif // CCfg_View_InheritanceIterator_h_INCLUDED
