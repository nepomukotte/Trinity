// =============================================================================
// InheritersIterator.h
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

#ifndef CCfg_View_InheritersIterator_h_INCLUDED
#define CCfg_View_InheritersIterator_h_INCLUDED

#include "InheritanceIterator.h"

namespace CCfg
{
namespace View
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * Enumerates every objects inheriting from the specified base object.
 *
 * @author fchateau
 */
class InheritersIterator : public InheritanceIterator
{
public:
	// construction/destruction
	InheritersIterator();
	explicit InheritersIterator(Attribute& baseObject);
	InheritersIterator(Attribute& baseObject, const Path& path);
	InheritersIterator(const InheritersIterator& value);
	virtual ~InheritersIterator();
	InheritersIterator& operator =(const InheritersIterator& value);

	// iteration
	bool hasNext();
};

}

}

#endif // CCfg_View_InheritersIterator_h_INCLUDED
