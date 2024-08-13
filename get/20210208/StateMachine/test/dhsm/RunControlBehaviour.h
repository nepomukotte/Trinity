// =============================================================================
// RunControlBehaviour.h
// -----------------------------------------------------------------------------
// enum dhsm::Behaviour
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_RunControlBehaviour_h_INCLUDED
#define dhsm_RunControlBehaviour_h_INCLUDED

namespace dhsm
{
enum Behaviour {
	Working,
	WillFail,
	WillCrash,
	WillHang,
	WillThrowCritical
};

} // namespace dhsm

#endif // dhsm_RunControlBehaviour_h_INCLUDED
