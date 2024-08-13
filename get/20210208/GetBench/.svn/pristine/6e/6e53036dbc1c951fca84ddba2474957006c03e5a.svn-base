/* =====================================================================================================================
 * dataRouter.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * dataRouter functions
 * Created on: Jun 20, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
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
#include "mdaq/DefaultPortNums.h"
#include "mdaq/utl/CmdLineArgs.h"
#include "mdaq/utl/Server.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"
#include "get/daq/DataRouter.h"
#include "mfm/FrameDictionary.h"
#include "utl/Logging.h"

#include <boost/algorithm/string.hpp>
namespace ba = boost::algorithm;
#include <memory>
#include <cstdlib>
#include <iostream>
#include <sstream>

using ::utl::net::SocketAddress;
using ::mdaq::utl::CmdLineArgs;
using ::mdaq::utl::Server;
using ::get::daq::DataRouter;

int main(int argc, char* argv[])
{
	// Setup logging backend
	::utl::BackendLogger::setBackend(::utl::LoggingBackendPtr(new ::mdaq::utl::ConsoleLoggingBackend));

	std::ostringstream usage;
	usage << "Usage: " << argv[0] << " <[ctrl-ipAddr]:[port]>  <[data-ipAddr]:[port]> [ICE|TCP|FDT] [ByteCounter|ByteStorage|FrameCounter|FrameStorage] [Ice arguments]";
	try
	{
		CmdLineArgs args = CmdLineArgs::build(argc, argv);

		if (args.size() > 1 and (args [1] == "-h" or args[1] == "--help")) throw usage.str();

		// Setting control endpoint
		SocketAddress ctrlEndpoint;
		if (args.size() > 1) ctrlEndpoint.fromString(args[1]);
		if (ctrlEndpoint.port() == 0) ctrlEndpoint.port() = ::mdaq::Default::dataRouterCtrlPortNum;
		if (not ctrlEndpoint.isValid()) ctrlEndpoint.ipAddress().in_addr() = INADDR_ANY;

		// Setting data flow endpoint (ip address will be same as control if not defined)
		SocketAddress flowEndpoint;
		if (args.size() > 2) flowEndpoint.fromString(args[2]);
		if (flowEndpoint.port() == 0) flowEndpoint.port() = ::mdaq::Default::dataRouterFlowPortNum;
		if (not flowEndpoint.isValid()) flowEndpoint.ipAddress().in_addr() = ctrlEndpoint.ipAddress().in_addr();

		// Setting data flow type
		std::string flowType = args.size() > 3 ? args[3] : "TCP";
		ba::to_upper(flowType);

		// Setting data processor type
		const std::string processorType = args.size() > 4 ? args[4] : "FrameStorage";

		// Load frame formats
		if ("FrameStorage" == processorType)
		{
			mfm::FrameDictionary::instance().addFormats("CoboFormats.xcfg");
		}

		// Disable IPv6 support
		// With Ice 3.5, support for IPv6 is enabled by default.
		// On sedipcg212 (Ubuntu 14.04), it causes "Address family is not supported by protocol" socket exceptions when using wildcard address.
		args.push_back("--Ice.IPv6=0");

		// Creating data router servant
		Server& server = Server::create(ctrlEndpoint, args);
		server.ic(); // Hack to ensure ICE communicator is initialized with desired arguments
		IceUtil::Handle< DataRouter > dataRouter(new DataRouter(flowEndpoint, flowType, processorType));
		server.addServant("DataRouter", dataRouter).start();
		dataRouter->runStart();
		server.waitForStop();
		return EXIT_SUCCESS;
	}
	catch (const std::exception & e) { LOG_FATAL() << e.what(); }
	catch (const std::string &  msg) { LOG_FATAL() << msg; }
	catch (const char*          msg) { LOG_FATAL() <<  msg; }
	catch (...)                      { LOG_FATAL() << "Unknown exception."; }
	return EXIT_FAILURE;
}
