// =============================================================================
// ServantApplication.h
// -----------------------------------------------------------------------------
// class dhsm::ServantApplication
// -----------------------------------------------------------------------------
// Frederic Chateau, Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay
// frederic.chateau@cea.fr, Shebli.Anvar@cea.fr
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef dhsm_ServantApplication_h_INCLUDED
#define dhsm_ServantApplication_h_INCLUDED

#include "Engine.h"
#include ICE_HEADER(ObjectAdapter.h)
#include ICE_HEADER(Communicator.h)
#include ICE_HEADER(Initialize.h)
#include <string>

namespace dhsm
{
/**
 * This class is the framework of the servant application of a distributed state
 * machine.
 */
class ServantApplication
{
public:
	ServantApplication();
	~ServantApplication();

	void init(int& argc, char* argv[], const std::string& endpoint);
	void init(int& argc, char* argv[], unsigned short port, const char* host = 0, int timeout = -1);
	void createServant(const char* identity, EnginePtr& engine);
	void run();

	Ice::ObjectAdapterPtr& objectAdapter() { return adapter; }
	Ice::CommunicatorPtr& communicator() { return ic; }

private:
	Ice::CommunicatorPtr ic;	///< the ICE communicator needed to use ICE.
	Ice::ObjectAdapterPtr adapter;	///< the object adapter where servants are registered.
};

}

#endif // dhsm_ServantApplication_h_INCLUDED
