/**
 * @file FrameSender.cpp
 * @date 31 oct. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetBench software project.
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


#include "mdaq/daq/TcpDataSender.hpp"
#include "mdaq/daq/IceDataSender.h"
#include "mdaq/utl/CmdLineArgs.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"
#include "mdaq/mt/Timer.h"
#include "mdaq/DefaultPortNums.h"
#include "mfm/FrameFormat.h"
#include "mfm/Frame.h"
#include "mfm/BitField.h"
#include "mfm/FrameDictionary.h"
#include "utl/Logging.h"

#ifdef VXWORKS
#include "mdaq/utl/net/ZBufConnection.h"
#include <IceE/IceE.h>
#include <IceE/Communicator.h>
#include <IceE/Time.h>

#else
#include <Ice/Ice.h>
#include <Ice/Communicator.h>
#include <IceUtil/Time.h>
#endif

#include <cmath>
#include <vector>
#include <fstream>

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
int frameSender(const char* argvargs)
{
	std::ostringstream usage;
	usage << "Usage: frameSender(\"ipAddr:[port] Ice|Tcp|ZBuf [test-duration-s]\")";
	CmdLineArgs args = CmdLineArgs::build("frameSender", argvargs);
#else
int main(int argc, char*argv[])
{
	std::ostringstream usage;
	usage << "Usage: " << argv[0] << " <[ipAddr]:[port]> <Ice|Tcp> [test-duration-s]";
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

		// Preparing data
		//const double PI = std::atan2(0.0f,-1.0f);
		std::auto_ptr<mfm::Frame> framePtr;
		try
		{
			LOG_INFO() << "Creating frame to send...";
			mfm::FrameDictionary::instance().addFormats("CoboFormats.xcfg");
			mfm::FrameFormat format = mfm::FrameDictionary::instance().findFormat(1, 5);
			std::auto_ptr<mfm::Frame> newFrame(format.createFrame()); // workaround for VxWorks
			framePtr = newFrame;
			framePtr->addItems(4*68*512);
			LOG_INFO() << "Initializing frame contents...";
			//size_t itemCount=0;
			mfm::Field field;
			mfm::BitField agetIdxField, chanIdxField, buckIdxField, sampleValueField;
			/*for (uint32_t agetId=0; agetId < 4; agetId++)
			{
				for (uint32_t timeId=0; timeId < 512; timeId++)
				{
					uint32_t value = (uint32_t) std::floor(500.+2500.*std::fabs(std::sin(2.*PI*timeId/512.)));
					for (uint32_t chanId = 0; chanId < 68; chanId++)
					{
						mfm::Item item = framePtr->itemAt(itemCount);
						if (itemCount <= 0)
						{
							field = item.field("");
							agetIdxField = field.bitField("agetIdx");
							chanIdxField = field.bitField("chanIdx");
							buckIdxField = field.bitField("buckIdx");
							sampleValueField = field.bitField("sample");
						}
						else
						{
							field = item.field(field);
							agetIdxField = field.bitField(agetIdxField);
							chanIdxField = field.bitField(chanIdxField);
							buckIdxField = field.bitField(buckIdxField);
							sampleValueField = field.bitField(sampleValueField);
						}
						agetIdxField.setValue(agetId);
						chanIdxField.setValue(chanId);
						buckIdxField.setValue(timeId);
						sampleValueField.setValue(value);
						itemCount++;
					}
				}
			}*/
		}
		catch (mfm::Exception & e)
		{
			LOG_FATAL() << "Could not create frame to send:" << e.what();
			return EXIT_FAILURE;
		}

		// Save frame in file for checks
		/*LOG_INFO() << "Saving frame";
		std::ofstream file("testFrame.graw", std::ios::binary | std::ios::trunc);
		framePtr->write(file);
		file.close();*/

		// Creating data sender
		std::auto_ptr < DataSender > dataSender;
		if (flowType == "Tcp")
		{
			LOG_DEBUG() << "Allocating TcpDataSender < TcpConnection >";
			dataSender.reset(new TcpDataSender < TcpConnection > (flowEndpoint));
		}
#ifdef VXWORKS
		else if (flowType == "ZBuf")
		{
			LOG_DEBUG() << "Allocating TcpDataSender < ZBufConnection >";
			dataSender.reset(new TcpDataSender < ZBufConnection > (flowEndpoint));
		}
#endif
		else if (flowType == "Ice")
		{
			LOG_DEBUG() << "Initializing Ice communicator for client access";
			ic = Ice::initialize(args);
			dataSender.reset(new IceDataSender(flowEndpoint, ic));
		}
		else throw "Unknown flow type: " + flowType;

		// Setting test duration
		int testDuration_s = 10; // 10 seconds by default
		if (args.size() > 3) testDuration_s = boost::lexical_cast<int>(args[3]);
		LOG_INFO() << "Test duration: " << testDuration_s << " s";

		// Sending data
		const size_t dataSize_B = framePtr->header().frameSize_B();
		LOG_INFO() << "Sending data frames of " << dataSize_B << " B during " << testDuration_s << " seconds";
		const uint8_t* dataBegin = reinterpret_cast<const uint8_t*>(framePtr->data());
		const uint8_t* dataEnd = dataBegin + dataSize_B+1;

		// Create circular buffer
		uint64_t dataBeginAddr = (uint64_t) &dataBegin;
		uint64_t dataEndAddr = (uint64_t) &dataEnd;
		mdaq::daq::CircularBuffer circularBuffer(dataBeginAddr, dataEndAddr);

		uint32_t eventIdx = 0;
		size_t nFramesSent = 0u;
		EndTestTimer endTestTimer;
		endTestTimer.start(testDuration_s);
		do
		{
			circularBuffer.reset();
			circularBuffer.writePtr.incr(dataSize_B);
			framePtr->headerField("eventIdx").setValue(eventIdx++);
			dataSender->flushData(circularBuffer);
			nFramesSent++;

		} while (not endTestTimer.endOfTest());

		LOG_INFO() << "Total data sent: " << nFramesSent << " frames in " << testDuration_s << " seconds";
		double const rate_Mbs = (8.*dataSize_B*nFramesSent/1000000.)/testDuration_s;
		LOG_INFO() << "Performance = " << rate_Mbs << " Mb/s";
		return EXIT_SUCCESS;
	}
	catch (const std::exception& e) { LOG_FATAL() << e.what();}
	catch (const std::string&  msg) { LOG_FATAL() << msg;}
	catch (const char*         msg) { LOG_FATAL() << msg;}
	catch (...)                     { LOG_FATAL() << "Unknown exception!";}
	return EXIT_FAILURE;
}



