// =============================================================================
// Visitor.cpp
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

#include "Visitor.h"

#include "Scalar.h"
#include "Object.h"

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Io::XmlFormat class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Destructor.
 */
Visitor::~Visitor()
{

}

/**
 * @fn bool Visitor::onObjectEnter(const Object& objectCfg)
 * Callback function invoked when we enter a new object while visiting the
 * configuration.
 * @param objectCfg the object we enter.
 * @return Returns false to indicate that visitation should ignore nested
 * attributes, true otherwise.
 */

/**
 * @fn bool Visitor::onObjectExit(const Object& objectCfg)
 * Callback function invoked when we leave an object while visiting the
 * configuration.
 * @param objectCfg the object we leave.
 * @return Returns false to indicate that visitation should ignore next
 * siblings.
 */

/**
 * @fn bool Visitor::onScalar(const Scalar& scalar)
 * Callback function invoked each time a scalar is encountered while
 * visiting the configuration.
 * @param scalar the scalar found.
 * @return Returns false to indicate that visitation should ignore next
 * siblings.
 */

}

