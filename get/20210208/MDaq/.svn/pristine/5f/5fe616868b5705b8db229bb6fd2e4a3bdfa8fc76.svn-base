/* =====================================================================================================================
 * eccServer.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * eccServer functions
 * Created on: Jan 11, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of the Mordicus Real-Time Software
 * Development Framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
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
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#include "mdaq/EccImpl.h"
#include "mdaq/DefaultPortNums.h"
#include "mdaq/utl/Server.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"
#ifdef HAVE_ICESTORM
#include "mdaq/utl/IceStormLoggingBackend.h"
#endif
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <cstdlib>

using ::mdaq::utl::Server;
using ::mdaq::EccImpl;
using ::mdaq::utl::CmdLineArgs;

static Server* eccServer = 0;

template< typename T >
bool getOptionValue(CmdLineArgs & args, const std::string & optionName, T & optionValue)
{
	CmdLineArgs::iterator arg;
	for (arg = args.begin(); arg != args.end(); arg++)
	{
		if (boost::starts_with(*arg, optionName + "="))
		{
			std::istringstream iss(arg->substr(optionName.size()+1));
			iss >> optionValue;
			args.erase(std::distance(args.begin(), arg));
			return static_cast<bool>(iss);
		}
	}
	return false;
}

int eccServerStop();

int eccServerRunOnPortWithArgs(int portNum, const CmdLineArgs& args = CmdLineArgs())
{
	try
	{
		eccServerStop();
		eccServer = &Server::create((uint16_t) portNum, args).addServant("Ecc", new EccImpl).start();
		eccServer->ic()->waitForShutdown();
		return EXIT_SUCCESS;
	}
	catch (const std::exception& e) { std::cerr << "*** " << e.what()  << std::endl; }
	catch (const std::string&  msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (const char*         msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (...)                     { std::cerr << "*** " << "UNKNOWN" << std::endl; }
	return EXIT_FAILURE;
}

int eccServerRunOnPort(int portNum = ::mdaq::Default::eccServerCtrlPortNum)
{
	return eccServerRunOnPortWithArgs(portNum);
}

int eccServerRun()
{
	return eccServerRunOnPort();
}

int eccServerStop()
{
	try
	{
		if (eccServer != 0)
		{
			eccServer->stop();
			eccServer = 0;
		}
		return 0;
	}
	catch (const std::exception& e) { std::cerr << "*** " << e.what()  << std::endl; }
	catch (const std::string&  msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (const char*         msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (...)                     { std::cerr << "*** " << "UNKNOWN" << std::endl; }
	return EXIT_FAILURE;
}

int main(int argc, char* argv[])
{
	using ::mdaq::utl::net::SocketAddress;
	try
	{
		CmdLineArgs args = CmdLineArgs::build(argc, argv);
		// Parse arguments for options
		if (std::find(args.begin(), args.end(), "-h") != args.end() or std::find(args.begin(), args.end(), "--help") != args.end())
		{
			std::cout << "Usage: " << argv[0] << " [[ipAddr][:ipPort]] [Ice arguments]" << std::endl;
			return EXIT_SUCCESS;
		}

		// Disable IPv6 support
		// With Ice 3.5, support for IPv6 is enabled by default.
		// On sedipcc212 (Ubuntu 14.04), it causes "Address family is not supported by protocol" socket exceptions when using wildcard address.
		args.push_back("--Ice.IPv6=0");

		// Set ICE communicator initialization data
		::mdaq::utl::Server::setInitData(args);

		// Setup logging backend
#ifdef HAVE_ICESTORM
		utl::BackendLogger::setBackend(utl::LoggingBackendPtr(new mdaq::utl::IceStormLoggingBackend("MDaqEccLogger")));
#else
		utl::BackendLogger::setBackend(utl::LoggingBackendPtr(new mdaq::utl::ConsoleLoggingBackend()));
#endif

		// Create server
		if (args.size() > 1)
		{
			mdaq::utl::net::SocketAddress sockAddr(args[1]);
			if (sockAddr.isValid())
			{
				args.erase(1);
				if (sockAddr.port() == 0)
				{
					sockAddr.port() = ::mdaq::Default::eccServerCtrlPortNum;
				}
				eccServer = &Server::create(sockAddr).addServant("Ecc", new EccImpl());
			}
			else
			{
				std::cerr << "--- Invalid IP address / port: " << sockAddr << std::endl;
				return EXIT_FAILURE;
			}
		}
		else
		{
			return eccServerRun();
		}

		// Start server thread
		eccServer->start();
		// Wait for shutdown
		eccServer->ic()->waitForShutdown();

		return EXIT_SUCCESS;
	}
	catch (const std::exception& e) { std::cerr << "*** " << e.what()  << std::endl; }
	catch (const std::string&  msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (const char*         msg) { std::cerr << "*** " <<   msg     << std::endl; }
	catch (...)                     { std::cerr << "*** " << "UNKNOWN" << std::endl; }
	return EXIT_FAILURE;
}
