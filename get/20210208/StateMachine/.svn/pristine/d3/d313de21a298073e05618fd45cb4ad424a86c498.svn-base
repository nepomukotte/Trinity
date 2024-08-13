// =============================================================================
// ServantApplication.cpp
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

#include "ServantApplication.h"
#include "EngineServant.h"
#include ICE_HEADER(Properties.h)
#include <iostream>
#include <sstream>

namespace dhsm
{
// =============================================================================
// dhsm::ServantApplication class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
ServantApplication::ServantApplication()
{
}

/**
 * Destructor.
 */
ServantApplication::~ServantApplication()
{
	if(ic.get() != 0)
	{
		ic->destroy();
	}
}

/**
 * Initializes the application.
 * @param argc the argc from "main".
 * @param argv the argv from "main".
 * @param endpoint the endpoint for ICE servants.
 */
void ServantApplication::init(int& argc, char* argv[], const std::string& endpoint)
{
	ic = Ice::initialize(argc, argv);
	adapter = ic->createObjectAdapterWithEndpoints("EngineServantAdapter", endpoint);
}

/**
 * Initializes the application.
 * @param argc the argc from "main".
 * @param argv the argv from "main".
 * @param port the port where the server listens.
 * @param host the host name where the server listens.
 * @param timeout the timeout for connections and calls.
 */
void ServantApplication::init(int& argc, char* argv[], unsigned short port, const char* host, int timeout)
{
	std::ostringstream endpoint;
	endpoint << "default -p " << port;
	if(host != 0)
	{
		endpoint << " -h " << host;
	}
	if(timeout >= 0)
	{
		endpoint << " -t " << timeout;
	}
	init(argc, argv, endpoint.str());
}

/**
 * Creates a servant engine for a distributed state machine.
 * @param identity the name of the ICE object implementing the state machine servant.
 * @param engine the engine executed by the servant (this Engine must be bound
 * to a StateMachine beforehand).
 */
void ServantApplication::createServant(const char* identity, EnginePtr& engine)
{
    std::auto_ptr<EngineServant> engineServant(new EngineServant);
    engineServant->setEngine(engine);
    Ice::ObjectPtr object = engineServant.release();

    adapter->add(object, ic->stringToIdentity(identity));
}

/**
 * Runs the application: executes the commands sent by the state machine master.
 */
void ServantApplication::run()
{
	adapter->activate();
    ic->waitForShutdown();
}

}
