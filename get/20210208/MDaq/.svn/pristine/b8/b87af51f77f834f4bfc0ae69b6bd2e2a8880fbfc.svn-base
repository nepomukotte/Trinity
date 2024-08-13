/**
 * @file PerfSend.cpp
 * @created Sep 14, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::daq::PerfSend
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

#include "mdaq/daq/TcpDataSender.hpp"
#include "mdaq/daq/IceDataSender.h"
#include "mdaq/utl/CmdLineArgs.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"
#include "mdaq/mt/Timer.h"
#include "mdaq/DefaultPortNums.h"
#include "utl/Logging.h"
#ifdef VXWORKS
#include "mdaq/utl/net/ZBufConnection.h"
#endif

#include ICE_HEADER(ICE_EDITION.h)
#include ICE_HEADER(Communicator.h)
#include ICE_HEADER(Config.h)
#ifdef ICEE_INT_VERSION
#include ICE_HEADER(Time.h)
#else
#include <IceUtil/Time.h>
#endif

#include <vector>

using ::utl::net::SocketAddress;
using ::mdaq::utl::CmdLineArgs;
using ::mdaq::utl::ConsoleLoggingBackend;
using ::mdaq::utl::net::TcpConnection;
using ::mdaq::daq::DataSender;
using ::mdaq::daq::TcpDataSender;
using ::mdaq::daq::IceDataSender;
using ::mdaq::mt::Timer;

#ifdef VXWORKS
using ::mdaq::utl::net::ZBufConnection;
#endif

class EndTestTimer : protected mdaq::mt::Timer
{
public:
	EndTestTimer() {}
	bool endOfTest() const { return  endOfTest_; }
	void start(int testDuration_s) { endOfTest_ = false; schedule(testDuration_s * 1000, false); }

protected:
	virtual void handler();

private:
	volatile bool endOfTest_;
};

void EndTestTimer::handler()
{
	endOfTest_ = true;
	PROMPT_MSG("End of test")
}

static Ice::CommunicatorPtr ic;

#ifdef VXWORKS
int perfsend(const char* argvargs)
{
	std::ostringstream usage;
	usage << "Usage: perfsend(\"ipAddr:[port] Ice|Tcp|ZBuf [data-buffer-size-KB] [test-duration-s]\")";
	CmdLineArgs args = CmdLineArgs::build("perfsend", argvargs);
#else
int main(int argc, char*argv[])
{
	std::ostringstream usage;
	usage << "Usage: " << argv[0] << " <[ipAddr]:[port]> <Ice|Tcp> [data-buffer-size-KB] [test-duration-s]";
	CmdLineArgs args = CmdLineArgs::build(argc, argv);
#endif
	try
	{
		::utl::BackendLogger::setBackend(::utl::LoggingBackendPtr(new ConsoleLoggingBackend));

		if (args.size() < 3) throw usage.str();

		// Setting flow endpoint
		SocketAddress flowEndpoint(args[1]);
		if (flowEndpoint.port() == 0) flowEndpoint.port() = ::mdaq::Default::dataRouterFlowPortNum;
		if (not flowEndpoint.isValid()) flowEndpoint.ipAddress().fromOctets(127,0,0,1);

		// Setting data flow type
		const std::string flowType = args[2];

		// Creating data sender
		std::auto_ptr < DataSender > dataSender;
		if (flowType == "Tcp")
		{
			LOG_INFO() << "Allocating TcpDataSender < TcpConnection >";
			dataSender.reset(new TcpDataSender < TcpConnection > (flowEndpoint));
		}
#ifdef VXWORKS
		else if (flowType == "ZBuf")
		{
			LOG_INFO() << "Allocating TcpDataSender < ZBufConnection >";
			dataSender.reset(new TcpDataSender < ZBufConnection > (flowEndpoint));
		}
#endif
		else if (flowType == "Ice")
		{
			LOG_INFO() << "Initializing Ice communicator for client access";
			ic = Ice::initialize(args);
			dataSender.reset(new IceDataSender(flowEndpoint, ic));
		}
		else throw "Unknown flow type: " + flowType;

		// Prepping data
		size_t dataSize_KB = 248*10;
		if (args.size() > 3) dataSize_KB = boost::lexical_cast<int>(args[3]);
		const size_t dataSize = dataSize_KB * 1024;
		std::vector<char> dataVect(dataSize+1); // allocate data
		char* data = &dataVect[0];
		char* dataEnd = data + dataSize + 1;
		LOG_INFO() << "Initializing " << dataSize_KB << " kibibytes of data...";
		for (size_t i=0u; i < dataSize; ++i) data[i] = i;

		// Create circular buffer
		uint64_t bufferStartAddr = (uint64_t) &data;
		uint64_t bufferEndAddr = (uint64_t) &dataEnd;
		mdaq::daq::CircularBuffer circularBuffer(bufferStartAddr, bufferEndAddr);

		// Setting test duration
		int testDuration_s = 10; // 10 seconds by default
		if (args.size() > 4) testDuration_s = boost::lexical_cast<int>(args[4]);

		// Sending data
		LOG_INFO() << "Sending data for " << testDuration_s << " seconds";
		uint64_t nBytesSent = 0u;
		EndTestTimer endTestTimer;
		endTestTimer.start(testDuration_s);

		do
		{
			circularBuffer.reset();
			circularBuffer.writePtr.incr(dataSize);
			dataSender->sendData(circularBuffer);
			nBytesSent += dataSize;
		} while (not endTestTimer.endOfTest());

		uint64_t nBytes_MB = nBytesSent/(1000000);
		LOG_INFO() << "Total data sent: " << nBytes_MB << " MB in " << testDuration_s << " seconds";
		LOG_INFO() << "Performance = " << (nBytes_MB*8.)/testDuration_s << " Mb/s";
		return EXIT_SUCCESS;
	}
	catch (const std::exception& e) { LOG_FATAL() << e.what(); }
	catch (const std::string&  msg) { LOG_FATAL() << msg; }
	catch (const char*         msg) { LOG_FATAL() << msg; }
	catch (...)                     { LOG_FATAL() << "Unknown exception!"; }
	return EXIT_FAILURE;

}
