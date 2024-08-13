// =============================================================================
// RunServantImpl.cpp
// -----------------------------------------------------------------------------
// class dhsm::RunServantImpl
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "RunServantImpl.h"
#include ICE_HEADER(ObjectAdapter.h)
#include ICE_HEADER(Communicator.h)
#include <iostream>

namespace dhsm
{
// =============================================================================
// dhsm::RunServantImpl class methods body
// -----------------------------------------------------------------------------

void RunServantImpl::shutdown(const Ice::Current& current)
{
	std::cout << "Shutting down" << std::endl;
	current.adapter->getCommunicator()->shutdown();
}

} // namespace dhsm
