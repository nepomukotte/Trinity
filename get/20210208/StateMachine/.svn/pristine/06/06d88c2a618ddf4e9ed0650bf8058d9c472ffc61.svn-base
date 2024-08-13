// =============================================================================
// RunServantImpl.h
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

#ifndef dhsm_RunServantImpl_h_INCLUDED
#define dhsm_RunServantImpl_h_INCLUDED

#include "dhsm/RunServant.h"

namespace dhsm
{
/**
 * This class implements the administrative interface for run servants that allow
 * the main test program to shutdown auxiliary test programs.
 */
struct RunServantImpl : public RunServant
{
	void shutdown(const Ice::Current& current);
};

} // namespace dhsm

#endif // dhsm_RunServantImpl_h_INCLUDED
