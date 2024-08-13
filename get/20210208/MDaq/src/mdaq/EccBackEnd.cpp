/*
 * $Id: EccBackEnd.cpp 1948 2017-11-23 16:14:28Z psizun $
 * @file EccBackEnd.cpp
 * @created 20 juin 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: patrick.sizun@cea.fr
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of
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

#include "mdaq/EccBackEnd.h"
#include "mdaq/cfg/NodeAtTop.h"
#include "mdaq/DefaultPortNums.h"
#include "utl/net/SocketAddress.h"
#include "utl/Logging.h"
#include "mdaq/rc/HardwareNode.h"
#include "mdaq/utl/Server.h"
#include "CCfg/Io/Url.h"

#include ICE_HEADER(Initialize.h)
#include ICE_HEADER(LocalException.h)
#include ICE_HEADER(Config.h)
#ifdef ICEE_INT_VERSION
#include ICE_HEADER(Thread.h)
#include ICE_HEADER(Time.h)
#else
#include <IceUtil/Thread.h>
#include <IceUtil/Time.h>
#endif

#include <boost/lexical_cast.hpp>
#include <stdexcept>
using ::utl::net::SocketAddress;

namespace mdaq {
//_________________________________________________________________________________________________
EccBackEnd::EccBackEnd() : currentHwNodeId_("")
{
	// Set options to communicate using Ice encoding 1.0
#if ICE_INT_VERSION >= 30500
	proxyOptions = "-e 1.0";
#endif
}
//_________________________________________________________________________________________________
void EccBackEnd::testConnectionToEcc()
{
	std::cout << "Ok" << std::endl;
}
//_________________________________________________________________________________________________
void EccBackEnd::testConnectionToHardware()
{
	if (hwNodes_.empty())
		throw utl::CmdException("No HW node has been defined.");

	for (HardwareNodeMap::const_iterator hwNodeIt = hwNodes_.begin(); hwNodeIt != hwNodes_.end(); ++hwNodeIt)
	{
		MSG() << "Testing connection to HW node " << hwNodeIt->first;
		hwNodeIt->second->hwServerCtrlPrx->testConnectionToHardware();
	}
}
//_________________________________________________________________________________________________
void EccBackEnd::shutdownHardwareServer()
{
	LOG_INFO() << "Shutting down hardware server...";
	currentHwNodePrx()->shutdown();
}
//_________________________________________________________________________________________________
void EccBackEnd::rebootHardware()
{
	LOG_INFO() << "Rebooting hardware...";
	hw::NodePrx proxyWithTimeout = hw::NodePrx::uncheckedCast(currentHwNodePrx()->ice_timeout(5000)); // ms
	try
	{
		proxyWithTimeout->reboot();
	}
	catch(const Ice::TimeoutException &)
	{
	}
}
//_________________________________________________________________________________________________
/**
 * Loads hardware description into current node.
 */
void EccBackEnd::loadHwDescription(const CCfg::CConfig & nodeCfg)
{
	cfg::NodeAtTop nodeDescr;
	try
	{
		// Parse config into a node description
		nodeDescr.fromCCfg(nodeCfg);
	}
	catch (CCfg::Io::Document::IoError& e) { LOG_ERROR() << e.what() ; throw utl::CmdException(e.what()) ; }
	catch (CCfg::CConfig::NotFound& e)     { LOG_ERROR() << e.what() ; throw utl::CmdException(e.what()) ; }
	catch (std::exception& e)              { LOG_ERROR() << e.what() ; throw utl::CmdException(e.what()) ; }
	catch (const char *msg)                { LOG_ERROR() << msg      ; throw utl::CmdException(msg)      ; }
	catch (const std::string& msg)         { LOG_ERROR() << msg      ; throw utl::CmdException(msg)      ; }
	catch (...)                            { LOG_ERROR() << "Unknown"; throw utl::CmdException("Unknown"); }

	// Loads the description to the target and instantiates the node devices
	// It first destroys the configuration loaded on the target
	currentHwNodePrx()->destroy();
	nodeDescr.id = std::string("HardwareNode@") + hwNodes_[currentHwNodeId_].endpoint.toString();
	currentHwNodePrx()->create(nodeDescr);

	// Update device structure of newly configured node.
	currentHwNode().updateDevices();
	LOG_DEBUG() << "Hardware description loaded.";
}
//_________________________________________________________________________________________________
/**
 * Loads hardware description into current node.
 */
void EccBackEnd::loadHwDescription(const std::string & cfgFile)
{
	LOG_INFO() << "Loading hardware node description from file '" << cfgFile << '\'';

	try
	{
		// Load contents of configuration file
		CCfg::Io::Document doc;
		doc.loadFromFile(cfgFile);
		CCfg::CConfig nodeCfg = doc.getConfig();

		// Load node description
		loadHwDescription(nodeCfg);
	}
	catch (utl::CmdException & e)          { throw; }
	catch (CCfg::Io::Document::IoError& e) { LOG_ERROR() << e.what(); throw utl::CmdException(e.what()) ; }
	catch (CCfg::CConfig::NotFound& e)     { LOG_ERROR() << e.what(); throw utl::CmdException(e.what()) ; }
	catch (std::exception& e)              { LOG_ERROR() << e.what(); throw utl::CmdException(e.what()) ; }
	catch (const char *msg)                { LOG_ERROR() << msg; throw utl::CmdException(msg)      ; }
	catch (const std::string& msg)         { LOG_ERROR() << msg; throw utl::CmdException(msg)      ; }
	catch (...)                            { LOG_ERROR() << "Unknown exception."; throw utl::CmdException("Unknown exception."); }
}
//_________________________________________________________________________________________________
void EccBackEnd::writeField(const std::string & registerName, const std::string& fieldName, Ice::Long fieldValue)
{
	const std::string deviceName = currentHwDeviceName();
	// LOG_DEBUG() << deviceName << "->writeField(" << registerName << ", " << fieldName << ", " << std::hex << std::showbase << fieldValue << ")" << std::dec;
	currentHwDevicePrx()->writeField(registerName, fieldName, fieldValue);
	// LOG_DEBUG() << deviceName << "->writeField() done";
}
//_________________________________________________________________________________________________
Ice::Long EccBackEnd::readField(const std::string& registerName, const std::string& fieldName)
{
	const std::string deviceName = currentHwDeviceName();
	// LOG_DEBUG() << deviceName << "->readField(" << registerName << ", " << fieldName << ")";
	Ice::Long valField = currentHwDevicePrx()->readField(registerName, fieldName);
	// LOG_DEBUG() << deviceName << "->readField() : " << std::hex << valField << std::dec;
	return valField;
}
//_________________________________________________________________________________________________
void EccBackEnd::addRegister(const ::mdaq::hw::RegisterConfig & cfg)
{
	currentHwDevicePrx()->addRegister(cfg);
}
//_________________________________________________________________________________________________
void EccBackEnd::addRegisters(const ::mdaq::hw::RegisterConfigList & cfg)
{
	currentHwDevicePrx()->addRegisters(cfg);
}
//_________________________________________________________________________________________________
void EccBackEnd::writeRegister(const std::string& registerName, Ice::Long regValue)
{
	const std::string deviceName = currentHwDeviceName();
	// LOG_DEBUG() << deviceName << "->writeRegister(" << registerName << ", " << std::hex << std::showbase << regValue << ")" << std::dec;
	currentHwDevicePrx()->writeRegister(registerName, regValue);
	// LOG_DEBUG() << deviceName << "->writeRegister() done";
}
//_________________________________________________________________________________________________
Ice::Long EccBackEnd::readRegister(const std::string& registerName)
{
	const std::string deviceName = currentHwDeviceName();
	// LOG_DEBUG() << deviceName << "->readRegister(" << registerName << ")";
	Ice::Long valReg = currentHwDevicePrx()->readRegister(registerName);
	// LOG_DEBUG() << deviceName << "->readRegister() : " << std::hex << valReg << std::dec;
	return valReg;
}
//_________________________________________________________________________________________________
/*
  * Sets data sender policy.
  * @param enable Whether to sent all data immediately or to optimize data chunk size.
  * @throws ::mdaq::utl::CmdException Raised no data sender currently exists.
  */
void EccBackEnd::setAlwaysFlushData(bool enable)
{
	LOG_DEBUG() << "Setting data sender 'alwaysFlushData' flag to " << enable;
	currentDaqCtrl()->setAlwaysFlushData(enable);
}
//_________________________________________________________________________________________________
/*
  * Selects which of the circular buffers should be read.
  * @param mask Mask of enabled buffers (1 bit per buffer).
  * @throws ::mdaq::utl::CmdException Raised no data sender currently exists.
  */
void EccBackEnd::setCircularBuffersEnabled(int mask)
{
	LOG_DEBUG() << "Enabling circular buffers with mask " << std::hex << std::showbase << mask << std::dec;
	currentDaqCtrl()->setCircularBuffersEnabled(mask);
}
//_________________________________________________________________________________________________
void EccBackEnd::daqStart()
{
	LOG_DEBUG() << "Calling daqStart";
	currentDaqCtrl()->daqStart();
	LOG_DEBUG() << "daqStart done.";
}
//_________________________________________________________________________________________________
void EccBackEnd::daqStop()
{
	LOG_DEBUG() << "Calling daqStop";
	LOG_INFO() << "daqStop";
	currentDaqCtrl()->daqStop();
	LOG_DEBUG() << "daqStop done.";
}
//_________________________________________________________________________________________________
void EccBackEnd::eccSleep(const size_t seconds)
{
	IceUtil::ThreadControl::sleep(IceUtil::Time::seconds(seconds));
}
//_________________________________________________________________________________________________
void EccBackEnd::connect(const std::string& objectType, const std::string & objectName, bool isBatched)
{
	try
	{
		if (objectType == "node")
		{
			// Connect to hardware node (and create a new HardwareNode instance).
			HardwareNodePtr hwNode = connectHwNode(objectName, isBatched);

			// Get string id of this new HW node
			const std::string hwNodeId = hwNode->hwServerCtrlPrx->name();
			// Add to known nodes
			addHwNode(hwNode);
			// Make new node the current node
			selectHwNodeById(hwNodeId);
		}
		else if (objectType == "device")
		{
			selectCurrentDevice(objectName);
		}
		else if (objectType == "daqserver" or objectType == "daq-ice")
		{
			connectRunProcessor(objectName, "ICE");
		}
		else
		{
			throw ::mdaq::utl::CmdException(std::string("Could not connect to '") + objectName +
					"': unknown type '" + objectType + '\'');
		}
	}
	catch (const ::mdaq::utl::CmdException & e)
	{
		throw;
	}
	catch (const Ice::ConnectionRefusedException & e)
	{
			throw ::mdaq::utl::CmdException(e.what());
	}
	catch (const Ice::Exception & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
	catch (const std::exception & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
}
//_________________________________________________________________________________________________
/**
 * Connects to an hardware node with given address and port.
 * @param endpoint IP address and port.
 */
HardwareNodePtr EccBackEnd::connectHwNode(const std::string & endpoint, bool isBatched)
{
	// Parse node endpoint into IP address and port
	utl::net::SocketAddress sockAddr(endpoint);
	if (0 == sockAddr.port())
		sockAddr.port() = Default::hwServerCtrlPortNum;

	LOG_INFO() << "Connecting to HW node at " << sockAddr << " with batch processing " << (isBatched?"enabled":"disabled");
	std::string hwNodeProxyString = sockAddr.iceProxyString("HwNode", "default", proxyOptions);

	Ice::ObjectPrx base = ::mdaq::utl::Server::ic()->stringToProxy(hwNodeProxyString);

	hw::NodePrx hwNodePrx = 0;
	try
	{
		hwNodePrx = hw::NodePrx::checkedCast(base);
	}
	catch (const Ice::SocketException & e)
	{
		throw ::mdaq::utl::CmdException(std::string("Connection to HW node \'") + sockAddr.toString() + "' failed:\n" + e.what());
	}
	if (not hwNodePrx)
	{
		LOG_ERROR() << "Connection to HW node '" << hwNodeProxyString << "' failed!";
		throw utl::CmdException(std::string("Connection to HW node '") + hwNodeProxyString + "' failed!");
	}
	return HardwareNodePtr(new rc::HardwareNode(sockAddr, hwNodePrx, isBatched));
}
//_________________________________________________________________________________________________
/**
 * Connects to DAQ control node for current hardware node and stores proxy.
 */
void EccBackEnd::connectDaqCtrl()
{
	currentHwNode().connectDaqCtrl();
}
//_________________________________________________________________________________________________
/**
 * Connects current node to its run processor.
 */
void EccBackEnd::connectRunProcessor(const std::string & dataRouterAddress, const std::string & dataRouterType)
{
	currentHwNode().connectRunProcessor(dataRouterAddress, dataRouterType);
}
//_________________________________________________________________________________________________
/**
 * Disconnects current node from its run processor.
 */
void EccBackEnd::disconnectRunProcessor()
{
	currentHwNode().disconnectRunProcessor();
}
//_________________________________________________________________________________________________
/**
 * Adds node to list of nodes known to ECC.
 */
void EccBackEnd::addHwNode(std::auto_ptr<rc::HardwareNode> hwNodePtr)
{
	const std::string hwNodeId = hwNodePtr->hwServerCtrlPrx->name();
	hwNodes_.insert(hwNodeId, hwNodePtr);
	hwNodes_[hwNodeId].updateDevices();

	LOG_DEBUG() << "Added HW node '" << hwNodes_[hwNodeId].hwServerCtrlPrx->name() << '\'';
}
//_________________________________________________________________________________________________
/**
 * Adds node to list of hardware nodes known to ECC.
 * @param endpoint Node socket address.
 */
void EccBackEnd::addHwNode(const ::utl::net::SocketAddress & endpoint)
{
	addHwNode(connectHwNode(endpoint.toString()));
}
//_________________________________________________________________________________________________
/**
 * Renames current node.
 * @param newName New name to give to the current node.
 */
void EccBackEnd::renameHwNode(const std::string & newName)
{
	// Rename node
	currentHwNodePrx()->setName(newName);
	LOG_DEBUG() << "HW node at '" << currentHwNode().endpoint << "' renamed into '" << newName << '\'';

	// Remove old node from map
	HardwareNodeMap::auto_type current_hw_node = hwNodes_.release(hwNodes_.find(currentHwNodeId_));
	HardwareNodePtr current_hw_node_ptr(current_hw_node.release());

	// Add node to map with new name
	currentHwNodeId_ = newName;
	hwNodes_.insert(currentHwNodeId_, current_hw_node_ptr);
	hwNodes_[currentHwNodeId_].updateDevices();
}
//_________________________________________________________________________________________________
void EccBackEnd::removeHwNode(const std::string & hwNodeId)
{
	if (hwNodes_.find(hwNodeId) != hwNodes_.end())
	{
		LOG_INFO() << "Removing HW node '" << hwNodeId << '\'';
		hwNodes_.erase(hwNodeId);
		if (hwNodeId == currentHwNodeId_)
			currentHwNodeId_ = "";
	}
}
//_________________________________________________________________________________________________
void EccBackEnd::removeAllHwNodes()
{
	hwNodes_.clear();
	currentHwNodeId_ = "";
}
//_________________________________________________________________________________________________
void EccBackEnd::selectHwNodeByName(const std::string& hwNodeName)
{
	HardwareNodeMap::const_iterator hwNodeIt = hwNodes_.begin();
	for (hwNodeIt = hwNodes_.begin(); hwNodeIt != hwNodes_.end(); ++hwNodeIt)
	{
		if (hwNodeIt->second->hwServerCtrlPrx->name() == hwNodeName)
		{
			currentHwNodeId_ = hwNodeIt->first;
			return;
		}
	}
	throw utl::CmdException(std::string("Could not find hardware node named '") + hwNodeName + "'!");
}
//_________________________________________________________________________________________________
void EccBackEnd::selectHwNodeByEndpoint(const ::utl::net::SocketAddress & endpoint)
{
	HardwareNodeMap::const_iterator hwNodeIt;
	for (hwNodeIt = hwNodes_.begin(); hwNodeIt != hwNodes_.end(); ++hwNodeIt)
	{
		if (hwNodeIt->second->endpoint.toString() == endpoint.toString())
		{
			currentHwNodeId_ = hwNodeIt->first;
			return;
		}
	}
	throw utl::CmdException(std::string("Could not find hardware node with endpoint '") + endpoint.toString() + "'!");
}
//_________________________________________________________________________________________________
void EccBackEnd::selectHwNodeById(const std::string& hwNodeId)
{
	HardwareNodeMap::const_iterator hwNodeIt = hwNodes_.find(hwNodeId);
	if (hwNodeIt == hwNodes_.end())
		throw utl::CmdException(std::string("Could not find hardware node '") + hwNodeId + "'!");
	currentHwNodeId_ = hwNodeIt->first;
}
//_________________________________________________________________________________________________
/**
 * Creates a device and selects it as the current device.
 */
void EccBackEnd::addDevice(const mdaq::hw::DeviceDescription & deviceDescription)
{
	currentHwNodePrx()->addDevice(deviceDescription);
	currentHwNode().updateDevices();
	selectCurrentDevice(deviceDescription.name);
}
//_________________________________________________________________________________________________
/**
 * Creates a device with registers, check write and read to registers and destroys the device.
 * Meant for testing the Mutant BRAM zone.
 */
std::string EccBackEnd::checkDevice(const ::mdaq::hw::DeviceDescription & descr, Ice::Long writeValue, Ice::Long firstRegOffset, Ice::Long regCount)
{
	// Create and select device
	addDevice(descr);

	// Add registers
	::mdaq::hw::RegisterConfigList registerCfg;
	registerCfg.resize(regCount);
	for (Ice::Long i=0; i < regCount; ++i)
	{
		registerCfg[i].descr.name = descr.name + "_Reg" + boost::lexical_cast< std::string >(i);
		registerCfg[i].descr.offset = firstRegOffset + i*descr.registerWidth;
	}
	addRegisters(registerCfg);

	// Loop over registers
	std::cout << std::hex;
	bool deviceOK = true;
	Ice::Long readValue = writeValue;
	Ice::Long regId;
	for (regId=0; regId < regCount; ++regId)
	{
		// Write and read register
		std::cout << "Offset 0x" << std::hex << registerCfg[regId].descr.offset << ": wrote 0x" << writeValue;
		writeRegister(registerCfg[regId].descr.name, writeValue);

		readValue = readRegister(registerCfg[regId].descr.name);
		std::cout << " read 0x" << std::hex << readValue << std::endl;

		if (readValue != writeValue) { deviceOK = false; break; }
	}

	// Destroy device
	removeDevice(descr.name);

	// Format result message
	std::ostringstream oss;
	if (not deviceOK)
	{
		oss << "FAILED: wrote 0x" << std::hex << writeValue << " to address 0x" <<  descr.baseAddress + registerCfg[regId].descr.offset << ", read back 0x" << readValue;
	}
	else
	{
		oss << "SUCCESS: successfully wrote and read back 0x" << std::hex << writeValue << " from address 0x" << (descr.baseAddress + registerCfg[0].descr.offset)
				<< " to address 0x" << (descr.baseAddress + registerCfg[regCount-1].descr.offset);
	}

	return oss.str();
}
//_________________________________________________________________________________________________
void EccBackEnd::removeDevice(const std::string & deviceName)
{
	LOG_DEBUG() << "Removing device '" << deviceName << "' in node '" << currentHwNode().endpoint << '\'';
	currentHwNodePrx()->removeDevice(deviceName);
	currentHwNode().updateDevices();
}
//_________________________________________________________________________________________________
void EccBackEnd::removeAllDevices()
{
	LOG_DEBUG() << "Removing all devices in node '" << currentHwNode().endpoint << '\'';
	currentHwNodePrx()->removeAllDevices();
	currentHwNode().updateDevices();
}
//_________________________________________________________________________________________________
/**
 * Selects a device from current node.
 * @param deviceName Name of the device.
 */
void EccBackEnd::selectCurrentDevice(const std::string & deviceName)
{
	// LOG_DEBUG() << "Selecting device '" << deviceName << "' in node '" << currentHwNode().endpoint << '\'';
	currentHwNode().selectDevice(deviceName);
}
//_________________________________________________________________________________________________
/**
 * Returns current node.
 */
rc::HardwareNode & EccBackEnd::currentHwNode()
{
	try
	{
		return hwNodes_.at(currentHwNodeId_);
	}
	catch (const boost::bad_ptr_container_operation &)
	{
		throw mdaq::utl::CmdException("No node has been selected.");
	}
}
//_________________________________________________________________________________________________
/**
 * Returns current node.
 */
const rc::HardwareNode & EccBackEnd::currentHwNode() const
{
	try
	{
		return hwNodes_.at(currentHwNodeId_);
	}
	catch (const std::out_of_range & e)
	{
		throw mdaq::utl::CmdException("No node has been selected.");
	}
}
//_________________________________________________________________________________________________
/**
 * Returns proxy to current node.
 */
hw::NodePrx& EccBackEnd::currentHwNodePrx()
{
	try
	{
		return hwNodes_.at(currentHwNodeId_).hwServerCtrlPrx;
	}
	catch (const boost::bad_ptr_container_operation &)
	{
	}
	throw mdaq::utl::CmdException("No node has been selected.");
}
//_________________________________________________________________________________________________
/**
 * Returns proxy to current device in current node.
 */
hw::DevicePrx& EccBackEnd::currentHwDevicePrx()
{
	return currentHwNode().currentDevicePrx();
}
//_________________________________________________________________________________________________
/**
 * Returns name of current device in current node.
 */
std::string EccBackEnd::currentHwDeviceName()
{
	return currentHwNode().currentDeviceName();
}
//_________________________________________________________________________________________________
/**
 * Returns proxy for DAQ control for current hardware node.
 */
mdaq::DaqCtrlNodePrx& EccBackEnd::currentDaqCtrl()
{
	return currentHwNode().daqCtrl();
}
//_________________________________________________________________________________________________
void EccBackEnd::getListOfHwNodes(hw::NodeList & list)
{
	list.clear();
	for (HardwareNodeMap::iterator hwNodeIt = hwNodes_.begin(); hwNodeIt != hwNodes_.end(); ++hwNodeIt)
	{
		list.push_back(hwNodeIt->second->hwServerCtrlPrx->name());
	}
}
//_________________________________________________________________________________________________
void EccBackEnd::getListOfDevices(hw::DeviceList & devices)
{
	currentHwNodePrx()->getListOfDevices(devices);
}
//_________________________________________________________________________________________________
void EccBackEnd::getListOfRegisters(hw::RegisterList & registers)
{
	currentHwDevicePrx()->getListOfRegisters(registers);
}
//_________________________________________________________________________________________________
void EccBackEnd::getListOfFields(const std::string & registerName, hw::FieldList & fields)
{
	currentHwDevicePrx()->getListOfFields(registerName, fields);
}
//_________________________________________________________________________________________________
/**
 * Returns name of currently selected hardware node.
 */
void EccBackEnd::getNameOfCurrentHwNode(std::string & nodeName)
{
	nodeName = currentHwNodePrx()->name();
}
//_________________________________________________________________________________________________
/**
 * Returns endpoint of currently selected hardware node.
 */
void EccBackEnd::getEndpointOfCurrentHwNode(std::string & nodeEndpoint)
{
	nodeEndpoint = currentHwNode().endpoint.toString();
}
//_________________________________________________________________________________________________
/**
 * Returns name of current device of currently selected hardware node.
 */
void EccBackEnd::getNameOfCurrentDevice(std::string & deviceName)
{
	deviceName = currentHwDeviceName();
}
//_________________________________________________________________________________________________
/**
 * Returns description of current device in currently selected hardware node.
 */
void EccBackEnd::getDeviceDescription(hw::DeviceDescription & descr)
{
	currentHwDevicePrx()->getDeviceDescription(descr);
}
//_________________________________________________________________________________________________
void EccBackEnd::getRegisterDescription(const std::string & registerName, hw::RegisterDescription & descr)
{
	currentHwDevicePrx()->getRegisterDescription(registerName, descr);
}
//_________________________________________________________________________________________________
void EccBackEnd::getFieldDescription(const std::string & registerName, const std::string & fieldName, hw::FieldDescription & descr)
{
	currentHwDevicePrx()->getFieldDescription(registerName, fieldName, descr);
}
//_________________________________________________________________________________________________
} /* namespace mdaq */
