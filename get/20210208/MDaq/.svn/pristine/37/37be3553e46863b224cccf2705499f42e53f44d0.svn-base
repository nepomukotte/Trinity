/**
 * @file HardwareNode.cpp
 * @date 8 oct. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors:
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

#include "mdaq/rc/HardwareNode.h"
#include "utl/net/SocketAddress.h"
#include "mdaq/DefaultPortNums.h"
#include "utl/Logging.h"
#include ICE_HEADER(Communicator.h)

namespace mdaq
{
namespace rc
{

HardwareNode::HardwareNode(::utl::net::SocketAddress address, const hw::NodePrx& proxy, bool isBatched)
	: endpoint(address), hwServerCtrlPrx(proxy), dataSenderCtrlPrx(0), currentRegProc_(regProcMap.end()), isBatched(isBatched)
{
}

HardwareNode::operator bool() const
{
	return this->hwServerCtrlPrx;
}

void HardwareNode::updateDevices()
{
	const std::string currentDeviceName = currentRegProc_ == regProcMap.end() ? "" : currentRegProc_->first;

	hw::DeviceMap deviceMap;
	hwServerCtrlPrx->getMapOfDevices(deviceMap);
	for (hw::DeviceMap::iterator devIt = deviceMap.begin(); devIt != deviceMap.end(); ++devIt)
	{
		regProcMap.insert(devIt->first, std::auto_ptr<hw::RegisterProcessor>(new hw::RemoteRegisterProcessor(devIt->second, isBatched)));
	}

	currentRegProc_ = regProcMap.find(currentDeviceName);
	if (currentRegProc_ == regProcMap.end() and not regProcMap.empty())
	{
		currentRegProc_ = regProcMap.begin();
	}
}

hw::RegisterProcessor::Map::iterator HardwareNode::findRegProc(const std::string& deviceName)
{
	hw::RegisterProcessor::Map::iterator regProcIt = regProcMap.find(deviceName);
	if (regProcIt == regProcMap.end())
		throw utl::CmdException(std::string("No device named '") + deviceName + "' in HW node '" + hwServerCtrlPrx->name() + "'!");

	return regProcIt;
}

hw::RegisterProcessor::Map::const_iterator HardwareNode::findRegProc(const std::string& deviceName) const
{
	hw::RegisterProcessor::Map::const_iterator regProcIt = regProcMap.find(deviceName);
	if (regProcIt == regProcMap.end())
		throw utl::CmdException(std::string("No device named '") + deviceName + "' in HW node '" + hwServerCtrlPrx->name() + "'!");

	return regProcIt;
}

void HardwareNode::selectDevice(const std::string& deviceName)
{
	currentRegProc_ = findRegProc(deviceName);
}

hw::RemoteRegisterProcessor& HardwareNode::regProc(const std::string& deviceName)
{
	return *dynamic_cast<hw::RemoteRegisterProcessor*>(findRegProc(deviceName)->second);
}

hw::RemoteRegisterProcessor& HardwareNode::currentRegProc()
{
	if (currentRegProc_ == regProcMap.end())
		throw ::mdaq::utl::CmdException("No device has been selected.");

	return *dynamic_cast<hw::RemoteRegisterProcessor*>(currentRegProc_->second);
}

hw::DevicePrx& HardwareNode::currentDevicePrx()
{
	return currentRegProc().devicePrx();
}

std::string HardwareNode::currentDeviceName()
{
	if (currentRegProc_ == regProcMap.end())
		throw ::mdaq::utl::CmdException("No device has been selected.");

	return currentRegProc_->first;
}

hw::DevicePrx& HardwareNode::findDevice(const std::string& deviceName)
{
	return regProc(deviceName).devicePrx();
}

//_________________________________________________________________________________________________
/**
 * Connects to node DAQ control and stores proxy.
 */
void HardwareNode::connectDaqCtrl()
{
	const ::utl::net::SocketAddress dataSenderAddress(endpoint.ipAddress(), mdaq::Default::dataSenderCtrlPortNum);
	std::string proxyOptions;
	// Set options to communicate using Ice encoding 1.0
#if ICE_INT_VERSION >= 30500
	proxyOptions = "-e 1.0";
#endif
	const std::string proxyString = dataSenderAddress.iceProxyString("DaqCtrlNode", "default", proxyOptions);
	LOG_DEBUG() << "Connecting DAQ control '" << proxyString << "' on node '" << hwServerCtrlPrx->name() << "'";
	Ice::ObjectPrx base = hwServerCtrlPrx->ice_getCommunicator()->stringToProxy(proxyString);

	mdaq::DaqCtrlNodePrx dataSenderNode = 0;
	try
	{
		dataSenderNode = mdaq::DaqCtrlNodePrx::checkedCast(base);
	}
	catch (const Ice::Exception & e)
	{
		throw utl::CmdException(std::string("Could not connect to ") + proxyString + ": " + e.what());
	}
	if (not dataSenderNode)
		throw utl::CmdException(std::string("Invalid proxy for DAQ control node: ") + proxyString);

	dataSenderCtrlPrx = dataSenderNode;
}
//_________________________________________________________________________________________________
/**
 * Returns node DAQ control proxy.
 */
mdaq::DaqCtrlNodePrx& HardwareNode::daqCtrl()
{
	if (not dataSenderCtrlPrx)
		throw utl::CmdException(std::string("No valid DAQ control node for node '") + hwServerCtrlPrx->name() + "'!");
	return dataSenderCtrlPrx;
}
//_________________________________________________________________________________________________
/**
 * Connects node to run processor.
 */
void HardwareNode::connectRunProcessor(const std::string & dataRouterAddress, const std::string & dataRouterType)
{
	try
	{
		connectDaqCtrl();
		LOG_DEBUG() << "Connecting to data router of type '" << dataRouterType << "' at address " << dataRouterAddress << " ...";
		daqCtrl()->connect(dataRouterType, dataRouterAddress);
		LOG_DEBUG() << "Connected.";
	}
	catch (const utl::CmdException & e)
	{
		throw utl::CmdException(std::string("Failed to connect to data router of type '") + dataRouterType + "': " + e.reason);
	}
	catch (const std::exception & e)
	{
		throw utl::CmdException(std::string("Failed to connect to data router of type '") + dataRouterType + "': " + e.what());
	}
}
//_________________________________________________________________________________________________
void HardwareNode::disconnectRunProcessor()
{
	if (dataSenderCtrlPrx)
	{
		try
		{
			dataSenderCtrlPrx->disconnect();
		}
		catch (const std::exception & /* e */)
		{
		}
	}
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace mdaq */
