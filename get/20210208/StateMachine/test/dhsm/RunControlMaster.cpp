// =============================================================================
// RunControlMaster.cpp
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

#include "RunControlMaster.h"
#include "RunServant.h"
#include "RunControlStateMachine.h"
#include ICE_HEADER(Process.h)
#include ICE_HEADER(LocalException.h)
#include <stdexcept>
#include <typeinfo>

namespace dhsm
{
// =============================================================================
// dhsm::RunControlMaster class methods body
// -----------------------------------------------------------------------------

RunControlMaster::RunControlMaster() : engine(), servantCount(1), servantHostIp("")
{
}

RunControlMaster::~RunControlMaster()
{
}

void RunControlMaster::init(int& argc, char* argv[])
{
	engine.init(argc, argv);

	parseArguments(argc, argv);

	StateMachinePtr sm(new RunControlStateMachine(Working));
	engine.setStateMachine(sm);
}

void RunControlMaster::start()
{
	engine.start();
}

void RunControlMaster::turnOn()
{
	engine.step("EvTurnOn");
}

void RunControlMaster::turnOff()
{
	engine.step("EvTurnOff");
}

void RunControlMaster::startRun()
{
	engine.step("EvStart");
}

void RunControlMaster::stopRun()
{
	engine.step("EvStop");
}

void RunControlMaster::pauseRun()
{
	engine.step("EvPause");
}

void RunControlMaster::resumeRun()
{
	engine.step("EvResume");
}

void RunControlMaster::run()
{
	try
	{
		connectServants();

		start();
		turnOn();
		startRun();
		pauseRun();
		resumeRun();
		stopRun();
		startRun();
		pauseRun();
		stopRun();
	}
	catch(const Exception& e)
	{
		std::cout << e.what() << std::endl;
	}
	turnOff();
	shutdownServants();
}

void RunControlMaster::parseArguments(int argc, char* argv[])
{
	for(int i = 1; i < argc-1; ++i)
	{
		if(not strcmp(argv[i], "-s"))
		{
			std::istringstream valueStr(argv[i+1]);
			valueStr >> servantCount;
			if(valueStr.bad())
			{
				std::cerr << "bad value for -s" << std::endl;
			}
		}
		else if(not strcmp(argv[i], "-h"))
		{
			std::istringstream valueStr(argv[i+1]);
			valueStr >> servantHostIp;
			if(valueStr.bad())
			{
				std::cerr << "bad value for -h" << std::endl;
			}    
		}
	}
}

void RunControlMaster::connectServants()
{
	for(int i = 1; i <= servantCount; ++i)
	{
		engine.connectServant("RunControlServant", 10000 + i, servantHostIp.size()==0 ? 0 : servantHostIp.c_str(), 3000);
	}
}

void RunControlMaster::shutdownServants()
{
	for(int i = 1; i <= servantCount; ++i)
	{
		try
		{
			Ice::CommunicatorPtr& ic = engine.communicator();
			std::ostringstream endpoint;
			endpoint << "RunServant:default -p " << 10000 + i << " -t 11000";
			if (servantHostIp.size()!=0)
			  endpoint << " -h " << servantHostIp;
			Ice::ObjectPrx servantBase = ic->stringToProxy(endpoint.str());
			RunServantPrx servant = RunServantPrx::checkedCast(servantBase);
			if(not servant)
			{
				throw std::runtime_error("cannot get servant administration interface");
			}
			else
			{
				servant->shutdown();
			}
		}
		catch(const Ice::SocketException& e)	// socket exceptions are not fatal: the servant may have crashed, etc.
		{
			std::cout << "Cannot close servant #" << i << std::endl;
		}
	}
}

} // namespace dhsm

// =============================================================================
// main
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	try
	{
		dhsm::RunControlMaster app;
		app.init(argc, argv);
		app.run();
	}
	catch(const Ice::Exception& e)
	{
		std::cerr << e << std::endl;
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
	return 0;
}
