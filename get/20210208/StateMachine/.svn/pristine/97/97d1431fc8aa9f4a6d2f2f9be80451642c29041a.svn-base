// =============================================================================
// ServantDescriptor.cpp
// -----------------------------------------------------------------------------
// class dhsm::ServantDescriptor
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "ServantDescriptor.h"

namespace dhsm
{
// =============================================================================
// dhsm::ServantDescriptor class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
ServantDescriptor::ServantDescriptor()
	: _status(Working), _errorMessage()
{
}

/**
 * Destructor.
 */
ServantDescriptor::~ServantDescriptor()
{
}

// =============================================================================
// dhsm::ServantDescriptor::Identity class methods body
// -----------------------------------------------------------------------------

ServantDescriptor::Identity::~Identity()
{
}

std::ostream& operator<<(std::ostream& stream, const ServantDescriptor::Identity& servantId)
{
	servantId.write(stream);
	return stream;
}

}
