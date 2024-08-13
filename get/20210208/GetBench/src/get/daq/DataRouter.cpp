/* =====================================================================================================================
 * DataRouter.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class get::daq::DataRouter
 * Created on: Jun 20, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * (c) Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GetBench package implementing the GET-specific
 * software modules for the MDaq framework.
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
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#include "DataRouter.h"
#include "mdaq/daq/IceDataReceiver.h"
#include "mdaq/daq/TcpDataReceiver.h"
#include "get/daq/FdtDataReceiver.h"
#ifndef __APPLE__
#include "mdaq/daq/NetPerfProcessor.h"
#include "FrameCounter.h"
#endif
#include "FrameStorage.h"
#include "ByteStorage.h"
#include "utl/Logging.h"
#include <boost/algorithm/string.hpp>
namespace ba = boost::algorithm;

using ::utl::net::SocketAddress;
using mdaq::daq::DataProcessorCore;
using mdaq::daq::DataReceiver;
using mdaq::daq::TcpDataReceiver;
using mdaq::daq::IceDataReceiver;

namespace get
{
namespace daq
{

DataRouter::DataRouter(const ::utl::net::SocketAddress& flowEndpoint, const std::string& flowType, const std::string& dataProcessorType)
{
	// Allocating DataReceiver corresponding to data flow type (Ice, Tcp...)
	createDataReceiver(flowEndpoint, flowType);

	// Allocating DataProcessor as specific DataProcessorCore and attaching it to data receiver
	std::auto_ptr < DataProcessorCore > dataProcessor;
	LOG_INFO() << "Creating DataProcessorCore of type '" << dataProcessorType << '\'';
	if (dataProcessorType == "ByteStorage")
	{
			dataProcessor.reset(new ByteStorage());
			dataReceiver->set_dataProcessorCore(dataProcessor);
	}
#ifndef __APPLE__
	else if (dataProcessorType == "ByteCounter")
	{
		::mdaq::daq::NetPerfProcessor* counter = new ::mdaq::daq::NetPerfProcessor();
		dataProcessor.reset(counter);
		dataReceiver->set_dataProcessorCore(dataProcessor);
		counter->startTimer();
	}
	else if (dataProcessorType == "FrameCounter" or dataProcessorType == "Debug")
	{
		FrameCounter* counter = new FrameCounter();
		dataProcessor.reset(counter);
		dataReceiver->set_dataProcessorCore(dataProcessor);
		counter->startTimer();
	}
#endif
	else if (dataProcessorType == "FrameStorage" or dataProcessorType == "DirectStorage")
	{
		dataProcessor.reset(new FrameStorage());
		dataReceiver->set_dataProcessorCore(dataProcessor);
	}
	else
	{
		throw "Unknown data processor type: '" + dataProcessorType + "'";
	}
}

DataRouter::DataRouter(const ::utl::net::SocketAddress& flowEndpoint, const std::string& flowType, std::auto_ptr < ::mdaq::daq::DataProcessorCore > dataProcessor)
{
	// Allocating DataReceiver corresponding to data flow type (Ice, Tcp...)
	createDataReceiver(flowEndpoint, flowType);

	// Allocating DataProcessor as specific DataProcessorCore and attaching it to data receiver
	dataReceiver->set_dataProcessorCore(dataProcessor);
}

void DataRouter::createDataReceiver(const ::utl::net::SocketAddress& flowEndpoint, const std::string& flowType)
{
	std::string dataFlowType = flowType;
	ba::to_upper(dataFlowType);
	LOG_INFO() << "Creating DataReceiver of type '" << flowType << "' with endpoint '" << flowEndpoint << '\'';
	// Allocating DataReceiver corresponding to data flow type (Ice, Tcp...)
	if (dataFlowType == "ICE")
	{
		dataReceiver.reset(new IceDataReceiver(flowEndpoint));
	}
	else if (dataFlowType == "TCP" or dataFlowType == "ZBUF")
	{
		dataReceiver.reset(new TcpDataReceiver(flowEndpoint, 0x100000));
	}
	else if (dataFlowType == "FDT")
	{
		dataReceiver.reset(new FdtDataReceiver(flowEndpoint, 0x100000));
	}
	else
	{
		throw "Unknown flow type: '" + dataFlowType + "'";
	}
}

DataRouter::~DataRouter()
{
}

void DataRouter::runConfig(const std::string&, const ::Ice::Current&)
{
	LOG_INFO() << "Configuring run processor...";
}

void DataRouter::runStart(const ::Ice::Current&)
{
	LOG_INFO() << "Starting run processor...";
	if (dataReceiver.get())
		dataReceiver->start();
}

void DataRouter::runStop(const ::Ice::Current&)
{
	LOG_INFO() << "Stopping run processor...";
	if (dataReceiver.get())
		dataReceiver->stop();
}

} // namespace daq
} // namespace get
