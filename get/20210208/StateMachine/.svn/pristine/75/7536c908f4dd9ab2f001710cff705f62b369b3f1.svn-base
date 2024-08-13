// =============================================================================
// RunControlServant.cpp
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

#include "RunControlServant.h"
#include "RunControlStateMachine.h"
#include "RunServantImpl.h"
#include "dhsm/Engine.h"
#include ICE_HEADER(ObjectAdapter.h)
#include ICE_HEADER(Communicator.h)
#include <sstream>
#include <iostream>

#ifdef RTEMS
#include <bsp.h>
#include <rtems/rtems_bsdnet.h>
#include "rtems_networkconfig.h"
#endif // RTEMS

namespace dhsm
{
// =============================================================================
// dhsm::RunControlServant class methods body
// -----------------------------------------------------------------------------

RunControlServant::RunControlServant()
	: app(), servantNumber(0), behaviour(Working), servantHostIp("")
{
}

RunControlServant::~RunControlServant()
{
}

void RunControlServant::run(int& argc, char* argv[])
{
	parseArguments(argc, argv);

	EnginePtr engine(new Engine);
	StateMachinePtr sm(new RunControlStateMachine(behaviour));
	engine->setStateMachine(sm);

	app.init(argc, argv, 10000 + servantNumber, servantHostIp.size()==0 ? 0 : servantHostIp.c_str());
	app.createServant("RunControlServant", engine);

	Ice::ObjectAdapterPtr& adapter = app.objectAdapter();
	Ice::CommunicatorPtr& ic = app.communicator();
	adapter->add(new RunServantImpl, ic->stringToIdentity("RunServant"));

	app.run();
}

void RunControlServant::parseArguments(int argc, char* argv[])
{
	for(int i = 1; i < argc-1; ++i)
	{
		if(not strcmp(argv[i], "-n"))
		{
			std::istringstream valueStr(argv[i+1]);
			valueStr >> servantNumber;
			if(valueStr.bad())
			{
				std::cerr << "bad value for -n" << std::endl;
			}
		}
		else if(not strcmp(argv[i], "-b") or not strcmp(argv[i], "--behaviour"))
		{
			if(not strcmp(argv[i+1], "work"))
			{
				behaviour = Working;
			}
			else if(not strcmp(argv[i+1], "fail"))
			{
				behaviour = WillFail;
			}
			else if(not strcmp(argv[i+1], "crash"))
			{
				behaviour = WillCrash;
			}
			else if(not strcmp(argv[i+1], "hang"))
			{
				behaviour = WillHang;
			}
			else if(not strcmp(argv[i+1], "critical"))
			{
				behaviour = WillThrowCritical;
			}
			else
			{
				std::cerr << "bad value for -b" << std::endl;
			}
		}
		else if(not strcmp(argv[i], "-h") or not strcmp(argv[i], "--host"))
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

} // namespace dhsm

// =============================================================================
// main
// -----------------------------------------------------------------------------

#ifdef RTEMS
int mainTask(int argc, char* argv[])
#else
int main(int argc, char* argv[])
#endif //RTEMS
{
	dhsm::RunControlServant app;
	app.run(argc, argv);
	return 0;
}

#ifdef RTEMS
/*
 * RTEMS Startup Task
 */
void *POSIX_Init(void *argument)
{
  char* servantNumber="1";
  char* behavior="work";
  char* ip="132.166.9.214";
  char* argv[] = {"RunControServant", "-n", servantNumber, "-b", behavior, "-h", ip};

  rtems_bsdnet_initialize_network ();

  mainTask(7,argv);
  exit (0);
}

#define CONFIGURE_APPLICATION_NEEDS_CONSOLE_DRIVER
#define CONFIGURE_APPLICATION_NEEDS_CLOCK_DRIVER

#define CONFIGURE_LIBIO_MAXIMUM_FILE_DESCRIPTORS 20
#define CONFIGURE_USE_IMFS_AS_BASE_FILESYSTEM

/* to be able to use posix mutex in RTEMS */
#define CONFIGURE_POSIX_INIT_THREAD_TABLE
#define CONFIGURE_MAXIMUM_POSIX_THREADS 100
#define CONFIGURE_MAXIMUM_POSIX_MUTEXES 100
#define CONFIGURE_MAXIMUM_POSIX_CONDITION_VARIABLES 300
#define CONFIGURE_MAXIMUM_POSIX_SEMAPHORES 300

#define CONFIGURE_EXECUTIVE_RAM_SIZE    (512*1024)
#define CONFIGURE_MAXIMUM_SEMAPHORES    50
#define CONFIGURE_MAXIMUM_TASKS         50

#define CONFIGURE_MICROSECONDS_PER_TICK 10000

#define CONFIGURE_INIT_TASK_STACK_SIZE  (10*1024)
#define CONFIGURE_INIT_TASK_PRIORITY    100
#define CONFIGURE_INIT_TASK_INITIAL_MODES (RTEMS_PREEMPT | \
                                           RTEMS_NO_TIMESLICE | \
                                           RTEMS_NO_ASR | \
                                           RTEMS_INTERRUPT_LEVEL(0))
#define CONFIGURE_INIT

#include <rtems/confdefs.h>

#endif // RTEMS

