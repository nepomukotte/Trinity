// =============================================================================
// Visitor.h
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

#ifndef CCfg_Visitor_h_INCLUDED
#define CCfg_Visitor_h_INCLUDED

namespace CCfg
{

//forward declarations
class Object;
class Scalar;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This is an interface to recursively visit each node of a configuration tree.
 *
 * @author fchateau
 */
class Visitor
{
public:
	// construction/destruction
	virtual ~Visitor();

	// visitation callbacks
	virtual bool onObjectEnter(const Object& objectCfg) = 0;
	virtual bool onObjectExit(const Object& objectCfg) = 0;
	virtual bool onScalar(const Scalar& scalar) = 0;
};

}

#endif // CCfg_Visitor_h_INCLUDED
