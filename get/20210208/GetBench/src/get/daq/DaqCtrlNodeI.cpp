/**
 * @file DaqCtrlNodeI.cpp
 * @date Sep 19, 2012
 * @author Shebli Anvar (shebli.anvar@cea.fr)
 * @author Joel Chavas (joel.chavas@cea.fr)
 *
 * @note class get::daq::DaqCtrlNodeI
 * @note Created on: Jun 21, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * @note This software is part of the GetBench package implementing the GET-specific software modules for the MDaq framework.
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
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "DaqCtrlNodeI.h"
#include "mdaq/daq/IceDataSender.h"
#include "mdaq/daq/TcpDataSender.hpp"
#include "get/daq/FakeDataSender.h"
#include "mdaq/daq/FdtDataSender.hpp"
#ifdef HAVE_MFM
#include "get/daq/FakeFrameSender.h"
#endif
#include "mdaq/mt/InterruptMonitorAdapter.hpp"
#include "mdaq/utl/CmdException.h"
#include "mdaq/utl/net/TcpConnection.h"
#include "mdaq/utl/net/ZBufConnection.h"
#include "mdaq/DefaultPortNums.h"
#include "utl/Logging.h"

#include <IceE/IceE.h>
#include <IceE/ObjectAdapter.h>

#include <logLib.h>
#include <sysLib.h>
#include <taskLib.h>
#include <types/vxTypesOld.h> // FUNCPTR
#include <hostLib.h> // hostGetByName

#include <boost/cstdint.hpp>

#include <unistd.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <exception>

using ::get::daq::MemRead;
using ::mdaq::mt::InterruptMonitorAdapter;
using ::mdaq::daq::DataSender;
using ::mdaq::daq::IceDataSender;
using ::mdaq::daq::StdDataSender;
using ::mdaq::daq::ZBufDataSender;
using ::mdaq::daq::FdtDataSender;
using ::mdaq::utl::net::TcpConnection;
using ::mdaq::utl::net::ZBufConnection;
using ::utl::net::SocketAddress;

namespace get
{
namespace daq
{

/**
 * Converts a time duration into a string using the most appropriate unit.
 * @param time Duration
 * @return String-ified version of the duration with unit.
 */
std::string stringify(const IceUtil::Time & time)
{
	std::ostringstream msg;
	IceUtil::Int64 seconds = time.toSeconds();
	if (seconds == 0) return "0 s";

	// Convert to days, hours, minutes, seconds
	IceUtil::Int64 days, hours, minutes;
	minutes = seconds/60;
	seconds = seconds%60;

	hours = minutes/60;
	minutes = minutes%60;

	days = hours/24;
	hours = hours%24;

	if (days > 0)
	{
		msg << days << " days";
	}
	if (hours > 0)
	{
		msg << ' ' << hours << " h";
	}
	if (minutes > 0)
	{
		msg << ' ' << minutes << " min";
	}
	if (seconds > 0)
	{
		msg << ' ' << seconds << " s";
	}
	return msg.str();
}
//_________________________________________________________________________________________________
const char* DaqCtrlNodeI::rateUnits[] = {" b/s", " kb/s", " Mb/s", " Gb/s"};
//_________________________________________________________________________________________________
DaqCtrlNodeI::DaqCtrlNodeI(cobo::CtrlNodeI & ctrlNode, int timerPeriod_s) : memReader(ctrlNode), timerPeriod_s(timerPeriod_s), isRunning_(false)
{
	createStatisticsTask();
}
//_________________________________________________________________________________________________
DaqCtrlNodeI::~DaqCtrlNodeI()
{
	deleteStatisticsTask();
}
//_________________________________________________________________________________________________
/**
 * Spawns a task to periodically print the data transfer statistics in the console.
 */
void DaqCtrlNodeI::createStatisticsTask()
{
	previousStats.reset();
	previousTime = IceUtil::Time::now();
	statsTaskId = taskSpawn("statsTask", 100, 0, 10000,
			reinterpret_cast<FUNCPTR>(statistics_printer),
			reinterpret_cast<int>(this), sysClkRateGet()*timerPeriod_s, 0, 0, 0, 0, 0, 0, 0, 0);
}
//_________________________________________________________________________________________________
void DaqCtrlNodeI::deleteStatisticsTask()
{
	taskDelete(statsTaskId);
}
//_________________________________________________________________________________________________
/**
 * Prints statistics periodically.
 * @param period_ticks Delay, in ticks, to wait before printing statistics again.
 */
void DaqCtrlNodeI::printStatistics(int period_ticks)
{
	if (UTL_LOG_LEVEL > 1) return;
	std::ostringstream msg;

	while (true)
	{
		// Wait ...
		taskDelay(period_ticks);
		// Protect task
		taskSafe();

		// Get new counts
		IceUtil::Time currentTime = IceUtil::Time::now();
		const MemRead::Stats currentStats = memReader.statistics();

		// Print rates
		if (currentStats.intrCount > previousStats.intrCount)
		{
			const double delay_s = (currentTime - previousTime).toSecondsDouble();
			double bitsPerSecond = double((currentStats.byteCount - previousStats.byteCount)*CHAR_BIT)/delay_s; // nb of bits per second
			// We do cannot compute the exact number of frames per second, just the number of interrupts that have been processed per second
			// const double framesPerSecond = double((currentStats.frameCount - previousStats.frameCount))/delay_s; // nb of frames per second
			// The total number of interrupts sent by the firmware gives a better estimate of the event rate
			const double intrPerSecond = (currentStats.intrCount - previousStats.intrCount)/delay_s;

			int iunit = 3;
			if      (bitsPerSecond >= 1e9)  { bitsPerSecond *= 1e-9;             }
			else if (bitsPerSecond >= 1e+6) { bitsPerSecond *= 1e-6; iunit--;    }
			else if (bitsPerSecond >= 1e+3) { bitsPerSecond *= 1e-3; iunit -= 2; }
			else                   {               iunit -= 3; }

			msg.str("");
			msg << ">>> " << std::fixed << std::setw(9) << std::setprecision(1) << bitsPerSecond << rateUnits[iunit] << ' ';
			if (intrPerSecond > 0) msg << std::fixed << std::setw(9) << std::setprecision(1) << intrPerSecond << " intr/s ";
			msg << std::setw(9) << currentStats.intrCount << " intr. ";
			msg << std::fixed << std::setw(6) << std::setprecision(2) << 100.*(1. - double(currentStats.memEmpty_B)/currentStats.memSize_B) << "% mem. occupied ";
			msg << " ; " << stringify(currentTime - runStartTime) << " since run start <<<";
			LOG_INFO() << msg.str();
		}

		// Update counts
		previousStats = currentStats;
		previousTime = currentTime;
	}
	// Release protection
	taskUnsafe();
}
//_________________________________________________________________________________________________
/**
 * Wrapper for non-static method print_statistics to be provided to taskSpawn.
 */
void DaqCtrlNodeI::statistics_printer(int obj, int arg)
{
	reinterpret_cast< DaqCtrlNodeI* >(obj)->printStatistics(arg);
}
//_________________________________________________________________________________________________
/**
 */
void DaqCtrlNodeI::dispatch(int imParam, int timeout_ms)
{
	mdaq::mt::InterruptMonitorAdapter<MemRead, MemRead>* im
		= reinterpret_cast<mdaq::mt::InterruptMonitorAdapter<MemRead, MemRead>*>(imParam);
	try
	{
		LOG_DEBUG() << "Starting interrupt monitoring with timeout of " << timeout_ms << " ms";
		im->start(timeout_ms);
	}
	catch (const ::mdaq::utl::CmdException& e)
	{
		LOG_ERROR() << e.reason;
	}
	catch (const std::exception& e)
	{
		LOG_ERROR() << e.what();
	}
	catch (const std::string & e)
	{
		LOG_ERROR() << e;
	}
	catch (const char* e)
	{
		LOG_ERROR() << e;
	}
	catch (...)
	{
		LOG_ERROR() << "Unknown exception when starting interrupt monitor";
	}
	im->stop();
	LOG_DEBUG() << "Monitoring of interrupts has been stopped";
}
//_________________________________________________________________________________________________
void DaqCtrlNodeI::createDataSender(const std::string & dataSenderType, const std::string & dataRouterRef, const Ice::Current & c)
{
	LOG_DEBUG() << "Creating data sender of type '" << dataSenderType << "'...";
	try
	{
		SocketAddress sockAddr(dataRouterRef);
		if (sockAddr.port() == 0) sockAddr.port() = mdaq::Default::dataRouterFlowPortNum;

		// Close previous connection
		memReader.resetDataSender();

		// Get address of local interface the TCP client should bind.
		SocketAddress bindAddress(INADDR_ANY, 0);
		const int dataIfAddr = hostGetByName("CoBoDataNode");
		if (ERROR != dataIfAddr)
		{
			bindAddress.ipAddress().in_addr() = dataIfAddr;
		}

		// Create new sender
		MemRead::DataSenderPtr dataSender;
		if (dataSenderType == "ICE") // Irfu Ice data connect
		{
			dataSender.reset(new IceDataSender(sockAddr, c.adapter->getCommunicator()));
		}
		else if (dataSenderType == "TCP") // Vanilla TCP connection
		{
			dataSender.reset(new StdDataSender(sockAddr, bindAddress));
		}
		else if (dataSenderType == "ZBUF") // Vanilla TCP connection from VxWorks ZBuf
		{
			dataSender.reset(new ZBufDataSender(sockAddr, bindAddress));
		}
		else if (dataSenderType == "FDT") // GANIL interface
		{
			dataSender.reset(new FdtDataSender<>(sockAddr, bindAddress));
		}
		else if (dataSenderType == "Debug") // Only count Bytes, no data sent
		{
			dataSender.reset(new FakeDataSender());
		}
#ifdef HAVE_MFM
		else if (dataSenderType == "DebugFrames") // Only count frames, no data sent
		{
			dataSender.reset(new FakeFrameSender());
		}
#endif
		else
		{
			std::string msg = std::string("Unknown data sender type '") + dataSenderType + '\'';
			LOG_ERROR() << msg;
			throw ::mdaq::utl::CmdException(msg);
		}
		memReader.resetDataSender(dataSender);

		LOG_INFO() << "Data sender of type '" << dataSenderType << "' successfully created.";
	}
	catch (const mdaq::utl::CmdException & e)
	{
		throw mdaq::utl::CmdException(e.reason);
	}
	catch (const Ice::Exception & e)
	{
		throw mdaq::utl::CmdException(e.what());
	}
	catch (const std::exception & e)
	{
		throw mdaq::utl::CmdException(e.what());
	}
}
//_________________________________________________________________________________________________
/**
 * Creates data sender and interrupt monitor.
 * @param dataRouterType Data flow type.
 * @param dataRouterRef Data router endpoint.
 * @param c Information about the current method invocation for servers.
 */
void DaqCtrlNodeI::connect(const std::string & dataRouterType, const std::string & dataRouterRef, const Ice::Current & c)
{
	try
	{
		// Create sender
		createDataSender(dataRouterType, dataRouterRef, c);

		// Create interrupt monitor
		LOG_DEBUG() << "Creating interrupt monitor";
		// InterruptMonitorAdapter(uint32_t baseAddress, int code,
		// 	C* actionObject, Action preWaitActionProc, Action postWaitActionProc, Action timeOutActionProc,
		//	P* handlerObject, Handler handlerProc)
		interruptMon.reset(new InterruptMonitorAdapter < MemRead, MemRead >
			(
				0x60040000u, 0,
				&memReader, 0, &MemRead::processInterrupt, &MemRead::sendHeartBeat,
				&memReader, &MemRead::onInterrupt
			)
		);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		throw mdaq::utl::CmdException(e.reason);
	}
	catch (const Ice::Exception & e)
	{
		throw mdaq::utl::CmdException(e.what());
	}
	catch (const std::exception & e)
	{
		throw mdaq::utl::CmdException(e.what());
	}
}
//_________________________________________________________________________________________________
void DaqCtrlNodeI::disconnect(const Ice::Current & /* c */)
{
	LOG_DEBUG() << "Disconnecting from data router";
	// Close previous connection
	memReader.resetDataSender();
}
//_________________________________________________________________________________________________
/*
 * Sets data sender policy.
 * @param enable Whether to sent all data immediately or to optimize data chunk size.
 * @throws ::mdaq::utl::CmdException Raised no data sender currently exists.
 */
void DaqCtrlNodeI::setAlwaysFlushData(bool enable, const Ice::Current &)
{
	memReader.setAlwaysFlushData(enable);
}
//_________________________________________________________________________________________________
/*
 * Selects which of the circular buffers should be read.
 * @param mask Mask of enabled buffers (1 bit per buffer).
 */
void DaqCtrlNodeI::setCircularBuffersEnabled(int asadMask, const Ice::Current &)
{
	memReader.setCircularBuffersEnabled(asadMask);
}
//_________________________________________________________________________________________________
/**
 * Starts data sender, updates circular buffer info and starts interrupt monitoring.
 */
void DaqCtrlNodeI::daqStart(const Ice::Current&)
{
	if (isRunning_) return;
	LOG_INFO() << "daqStart";
	// Start data sender
	try
	{
		memReader.dataSender().start();
		memReader.sendTopology();
	}
	catch (const std::exception & e)
	{
		throw ::mdaq::utl::CmdException(std::string("Could not establish data link. ") + e.what());
	}
	runStartTime = IceUtil::Time::now();
	// Reset circular buffer
	memReader.resetCircularBuffers();
	previousStats.reset();

	// Start interrupt monitoring thread
	const int heartBeatPeriod_ms = 3000;
	taskSpawn("intrTask", 110, 0, 20000, reinterpret_cast< FUNCPTR >(&dispatch),
			reinterpret_cast<int>(this->interruptMon.get()), heartBeatPeriod_ms,
			0,0,0,0,0,0,0,0);
	isRunning_ = true;
	// Wait
	IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(20)); // IceUtil::Time::milliSeconds(1000));
}
//_________________________________________________________________________________________________
/**
 * Stops monitoring of interrupts and sends all data left in buffer.
 */
void DaqCtrlNodeI::daqStop(const Ice::Current &)
{
	LOG_INFO() << "daqStop";
	isRunning_ = false;
	if (interruptMon.get())
	{
		// Stop monitoring interrupts
		interruptMon->stop();
		// Wait for last interrupt to be processed (e.g. for data to be sent)
		for (size_t numTries=0; numTries < 10 and not interruptMon->waitForStopped(1000); numTries++)
		{
			LOG_INFO() << "Waiting for interrupt monitor to stop...";
		}
	}
	flushData();
	memReader.dataSender().stop();
}
//_________________________________________________________________________________________________
/**
 * Flushes all data contained in circular buffer.
 */
void DaqCtrlNodeI::flushData(const Ice::Current &)
{
	LOG_INFO() << "Flushing data...";
	memReader.flushData();
}
//_________________________________________________________________________________________________
void DaqCtrlNodeI::shutdown(const Ice::Current& c)
{
	LOG_DEBUG() << "Shutting down DAQ control node";
    c.adapter->getCommunicator()->shutdown();
}
//_________________________________________________________________________________________________
} // namespace daq
} // namespace get
