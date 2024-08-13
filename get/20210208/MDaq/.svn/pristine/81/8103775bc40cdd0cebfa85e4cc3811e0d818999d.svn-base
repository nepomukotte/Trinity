/**
 * @file PerfRcve.cpp
 * @created Sep 14, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::daq::PerfRcve
 * Created on: Sep 14, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the MDaq software project.
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

#include "mdaq/daq/TcpDataReceiver.h"
#include "mdaq/daq/IceDataReceiver.h"
#include "mdaq/DefaultPortNums.h"
#include "mdaq/daq/NetPerfProcessor.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"

#include <boost/lexical_cast.hpp>
#include <memory>


using ::utl::net::SocketAddress;
using ::mdaq::utl::CmdLineArgs;
using ::mdaq::utl::ConsoleLoggingBackend;
using ::mdaq::daq::DataReceiver;
using ::mdaq::daq::DataProcessorCore;
using ::mdaq::daq::TcpDataReceiver;
using ::mdaq::daq::IceDataReceiver;
using ::mdaq::daq::NetPerfProcessor;

int main(int argc, char*argv[])
{
	std::ostringstream usage;
	usage << "Usage: " << argv[0] << " <[ipAddr]:[port]> <Ice|Tcp> [bufSize]";
	try
	{
		::utl::BackendLogger::setBackend(::utl::LoggingBackendPtr(new ConsoleLoggingBackend));
		CmdLineArgs args = CmdLineArgs::build(argc, argv);

		if (args.size() < 3) throw usage.str();

		// Setting flow endpoint
		SocketAddress flowEndpoint(args[1]);
		if (flowEndpoint.port() == 0) flowEndpoint.port() = ::mdaq::Default::dataRouterFlowPortNum;
		if (not flowEndpoint.isValid()) flowEndpoint.ipAddress().fromOctets(127, 0, 0, 1);

		// Setting data flow type
		const std::string flowType = args[2];

		// Creating data receiver
		std::auto_ptr < DataReceiver > dataReceiver;
		int bufSize = 0x100000;
		if (flowType == "Tcp")
		{
			if (args.size() > 3) bufSize = boost::lexical_cast<size_t>(args[3]);
			MSG() << "Alloc. TcpDataReceiver @ " << flowEndpoint.toString() << " with " << bufSize << " bytes buffer size";
			dataReceiver.reset(new TcpDataReceiver(flowEndpoint, bufSize));
		}
		else if (flowType == "Ice")
		{
			MSG() << "Allocating IceDataReceiver @ " << flowEndpoint.toString() << " with " << bufSize << " bytes buffer size";
			dataReceiver.reset(new IceDataReceiver(flowEndpoint));
		}
		else throw "Unknown flow type: " + flowType;

		// Setting DataProcessorCore
		NetPerfProcessor* netPerf = new NetPerfProcessor;
		std::auto_ptr < DataProcessorCore > netPerf_(netPerf);
		dataReceiver->set_dataProcessorCore(netPerf_);

		// Launching
		netPerf->startTimer();
		dataReceiver->start();
		dataReceiver->waitForStop();

		return EXIT_SUCCESS;
	}
	catch (const std::exception& e) { LOG_FATAL() << e.what() ; }
	catch (const std::string&  msg) { LOG_FATAL() << msg      ; }
	catch (const char*         msg) { LOG_FATAL() << msg      ; }
	catch (...)                     { LOG_FATAL() << "UNKNOWN"; }
	return EXIT_FAILURE;

}
