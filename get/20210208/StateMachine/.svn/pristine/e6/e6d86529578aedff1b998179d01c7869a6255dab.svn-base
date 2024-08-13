// =============================================================================
// RunControlMaster.h
// -----------------------------------------------------------------------------
// class dhsm::RunControlMaster
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_RunControlMaster_h_INCLUDED
#define dhsm_RunControlMaster_h_INCLUDED

#include "dhsm/MasterEngine.h"
#include ICE_HEADER(Communicator.h)
#include <string>

namespace dhsm
{
/**
 * This is a test class that uses the state-machine framework to implement a
 * simple run control.
 */
class RunControlMaster
{
public:
	RunControlMaster();
	~RunControlMaster();

	void init(int& argc, char* argv[]);
	void run();

private:
	void parseArguments(int argc, char* argv[]);
	void connectServants();
	void shutdownServants();
	void start();
	void turnOn();
	void turnOff();
	void startRun();
	void stopRun();
	void pauseRun();
	void resumeRun();

	MasterEngine engine;
	int servantCount;
	std::string servantHostIp;
};

} // namespace dhsm

#endif // dhsm_RunControlMaster_h_INCLUDED
