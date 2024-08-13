// =============================================================================
// RunControlServant.h
// -----------------------------------------------------------------------------
// class dhsm::RunControlServant
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_RunControlServant_h_INCLUDED
#define dhsm_RunControlServant_h_INCLUDED

#include "dhsm/ServantApplication.h"
#include "RunControlBehaviour.h"
#include <string>

namespace dhsm
{
/**
 * This is a test class simulating the state-machine of a data acquisition node
 * controlled by a run control.
 */
class RunControlServant
{
public:
	RunControlServant();
	~RunControlServant();

	void run(int& argc, char* argv[]);

private:
	void parseArguments(int argc, char* argv[]);

	ServantApplication app;
	int servantNumber;
	Behaviour behaviour;
	std::string servantHostIp;
};

} // namespace dhsm

#endif // dhsm_RunControlServant_h_INCLUDED
