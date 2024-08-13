/**
 * @file ThreadSafeLoggerTest.cpp
 * @date 12 May 2016
 * @author anvar
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun (patrick.sizun@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 *
 * This file is part of the MDaq software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "mdaq/utl/CmdLineArgs.h"

#ifdef VXWORKS

#include <vxWorks.h>
#include <taskLib.h>
#include <semLib.h>
#else
#include <boost/thread.hpp>
#endif
#include <boost/lexical_cast.hpp>
#include <utl/Logging.h>
#include <algorithm>
#include <iostream>
#include <vector>

using ::mdaq::utl::CmdLineArgs;
using ::utl::LoggingBackendPtr;
using ::utl::BackendLogger;
using ::utl::StdoutLoggingBackend;

#ifdef VXWORKS
void logMessages(SEM_ID joinSem)
#else
void logMessages()
#endif
{
	LOG_DEBUG() << "This is a DEBUG message.";
	LOG_INFO() << "This is a INFO message.";
	LOG_WARN() << "This is a WARN message.";
	LOG_ERROR() << "This is a ERROR message.";
	LOG_FATAL() << "This is a FATAL message.";

#ifdef VXWORKS
	semGive(joinSem);
#endif
}

#ifdef VXWORKS
int threadsafe_logger_test(const char* argvargs)
{
	std::ostringstream usage;
	usage << "Usage: threadsafe_logger_test";
	CmdLineArgs args = CmdLineArgs::build("threadsafe_logger_test", argvargs);
#else
int main(int argc, char*argv[])
{
	std::cout << std::string(80, '-') << std::endl;

	std::ostringstream usage;
	usage << "Usage: " << argv[0] << " [number of threads]";
	CmdLineArgs args = CmdLineArgs::build(argc, argv);
#endif

	try
	{
		size_t nThreads = 4;

		// Parse arguments for options
		if (std::find(args.begin(), args.end(), "-h") != args.end() or std::find(args.begin(), args.end(), "--help") != args.end())
		{
			std::cout << usage.str() << std::endl;
			return EXIT_SUCCESS;
		}
		else if (args.size() > 1)
		{
			nThreads = boost::lexical_cast<size_t>(args[1]);
		}

		// Setup logging backend
		BackendLogger::setBackend(LoggingBackendPtr(new StdoutLoggingBackend));

#ifdef VXWORKS

std::vector<int> threads(nThreads, -1);
		std::vector<SEM_ID> joinSems(nThreads, NULL);
		int lowerPrio = -1;
		taskPriorityGet(taskIdSelf(), &lowerPrio);
		lowerPrio += 1; // Threads are launched with lower priority than main
		for (size_t i=0u; i < nThreads; ++i)
		{
			joinSems[i] = semBCreate(SEM_Q_FIFO, SEM_EMPTY);
			threads[i] = taskSpawn(NULL, lowerPrio, 0, 10000, (FUNCPTR) logMessages, (int) joinSems[i],1,2,3,4,5,6,7,8,9);
		}
		for (size_t i=0u; i < nThreads; ++i)
		{
			semTake(joinSems[i], WAIT_FOREVER);
		}

#else
		boost::thread_group thg;
		for (size_t i=0u; i < nThreads; ++i)
		{
			thg.create_thread(logMessages);
		}
		thg.join_all();

		#endif

		return EXIT_SUCCESS;
	}
	catch (const std::exception& e) { LOG_FATAL() << e.what(); }
	catch (const std::string&  msg) { LOG_FATAL() << msg; }
	catch (const char*         msg) { LOG_FATAL() << msg; }
	catch (...)                     { LOG_FATAL() << "Unknown exception!"; }

	return EXIT_FAILURE;
}
