/**
 * @file IceStormLoggingBackendTest.cpp
 * @date 12 déc. 2012
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
#include "mdaq/utl/IceStormLoggingBackend.h"
#include <cstdlib>
#include <algorithm>
#include <iostream>

using ::mdaq::utl::CmdLineArgs;

#ifdef VXWORKS
int icestorm_logger_test_publisher(const char* argvargs)
{
	std::ostringstream usage;
	usage << "Usage: icestorm_logger_test_publisher(\"--Ice.Config=[ICE_CONFIG_FILE] [message] [FATAL|ERROR|WARN|INFO|DEBUG] [topic-name]\")";
	CmdLineArgs args = CmdLineArgs::build("icestorm_logger_test_publisher", argvargs);
#else
int main(int argc, char*argv[])
{
	std::ostringstream usage;
	usage << "Usage: " << argv[0] << "--Ice.Config=[ICE_CONFIG_FILE] [message] [FATAL|ERROR|WARN|INFO|DEBUG] [topic-name]";
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

		// Set ICE communicator initialization data
		std::cout << "Setting up ICE communicator initialization data" << std::endl;
		::mdaq::utl::Server::setInitData(args);

		// Parse options
		std::string msg = "This is a test";
		if (args.size() > 2)
			msg = args.at(2);
		std::string level = "INFO";
		if (args.size() > 3)
			level = args.at(3);
		std::string topic = "TestLogger";
		if (args.size() > 4)
			topic = args.at(4);

		std::cout << "Instantiating IceStorm logging backend" << std::endl;
		utl::BackendLogger::setBackend(utl::LoggingBackendPtr(new mdaq::utl::IceStormLoggingBackend(topic)));

		std::cout << "Sending message" << std::endl;
		if (level == "FATAL")
		{
			LOG_FATAL() << msg;
		}
		else if (level == "ERROR")
		{
			LOG_ERROR() << msg;
		}
		else if (level == "WARN")
		{
			LOG_WARN() << msg;
		}
		else if (level == "DEBUG")
		{
			LOG_DEBUG() << msg;
		}
		else
		{
			LOG_INFO() << msg;
		}

		return EXIT_SUCCESS;
	}
	catch (const Ice::Exception & e) { LOG_FATAL() << e.what(); }
	catch (const std::exception& e) { LOG_FATAL() << e.what(); }
	catch (const std::string&  msg) { LOG_FATAL() << msg; }
	catch (const char*         msg) { LOG_FATAL() << msg; }
	catch (...)                     { LOG_FATAL() << "Unknown exception!"; }

	return EXIT_FAILURE;
}



