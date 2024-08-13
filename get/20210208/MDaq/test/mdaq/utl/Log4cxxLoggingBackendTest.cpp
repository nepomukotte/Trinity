/**
 * @file Log4cxxLoggingBackendTest.cpp
 * @date 13 déc. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
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
#include "mdaq/utl/Server.h"
#include "mdaq/utl/Log4cxxLoggingBackend.h"
#include <log4cxx/basicconfigurator.h>
#include <algorithm>
#include <iostream>

using ::mdaq::utl::CmdLineArgs;

#ifdef VXWORKS
int log4cxx_logger_test(const char* argvargs)
{
	std::ostringstream usage;
	usage << "Usage: log4cxx_logger_test";
	CmdLineArgs args = CmdLineArgs::build("log4cxx_logger_test", argvargs);
#else
int main(int argc, char*argv[])
{
	std::ostringstream usage;
	usage << "Usage: " << argv[0];
	CmdLineArgs args = CmdLineArgs::build(argc, argv);
#endif

	try
	{
		// Parse arguments for options
		if (std::find(args.begin(), args.end(), "-h") != args.end() or std::find(args.begin(), args.end(), "--help") != args.end())
		{
			std::cout << usage.str() << std::endl;
			return EXIT_SUCCESS;
		}

		// Set up a simple configuration that logs on the console.
		log4cxx::BasicConfigurator::configure();

		// Setup logging backend
		utl::BackendLogger::setBackend(utl::LoggingBackendPtr(new mdaq::utl::Log4cxxLoggingBackend("MDaq")));

		LOG_DEBUG() << "This is a DEBUG message.";
		LOG_INFO() << "This is a INFO message.";
		LOG_WARN() << "This is a WARN message.";
		LOG_ERROR() << "This is a ERROR message.";
		LOG_FATAL() << "This is a FATAL message.";

		return EXIT_SUCCESS;
	}
	catch (const Ice::Exception & e) { LOG_FATAL() << e.what(); }
	catch (const std::exception& e) { LOG_FATAL() << e.what(); }
	catch (const std::string&  msg) { LOG_FATAL() << msg; }
	catch (const char*         msg) { LOG_FATAL() << msg; }
	catch (...)                     { LOG_FATAL() << "Unknown exception!"; }

	return EXIT_FAILURE;
}


