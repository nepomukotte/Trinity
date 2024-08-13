/* =====================================================================================================================
 * BackEnd.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class get::rc::BackEnd
 * Created on: Apr 4, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GET data acquisition framework.
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

#include "get/rc/BackEnd.h"
#include "get/mt/AlarmLogger.h"
#include "get/rc/ConfigManager.h"
#include "get/rc/ConfigFileManager.h"
#if CFG_DB_SUPPORT_ENABLED
#include "get/rc/ConfigDatabaseManager.h"
#endif
#include "get/rc/CoBoPLLManager.h"
#include "get/asad/utl.h"
#include "get/mutant/AlignmentServer.h"
#include "mdaq/DefaultPortNums.h"
#include "mdaq/utl/CmdException.h"
#include "mdaq/utl/Server.h"
#include "dhsm/State.h"
#include "dhsm/StateMachine.h"
#include "dhsm/Message.h"
#include <CCfg/Io/BufferStorage.h>
#include <CCfg/Io/Document.h>
#include <CCfg/Io/Url.h>
#include <CCfg/CConfig.h>
#include <CCfg/Object.h>
#include <CCfg/View/InheritersIterator.h>
#include <CCfg/Io/XmlFormat.h>
#include "CCfg/Io/IncludeHandler.h"
#include "utl/Logging.h"
#include ICE_HEADER(Initialize.h)
#include ICE_HEADER(LocalException.h)
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/classification.hpp>
namespace ba = boost::algorithm;
#include <sstream>
#include <iomanip>

using ::mdaq::utl::CmdException;
using ::mdaq::utl::Server;
using namespace dhsm;
using namespace CCfg;

#define CATCH_SM_EXCEPTIONS(state)						\
	catch(const CmdException & e)			\
	{													\
		LOG_ERROR() << e.reason;						\
		handleStateMachineException(state, e.reason);	\
	}													\
	catch (const Ice::LocalException & e) 				\
	{													\
		LOG_ERROR() << e.what();						\
		handleStateMachineException(state, std::string("Ice exception: ") + e.what());	\
	}													\
	catch(const std::exception & e)						\
	{													\
		LOG_ERROR() << e.what();						\
		handleStateMachineException(state, std::string("Exception: ") + e.what());	\
	}													\
	catch (const std::string & msg)						\
	{													\
		LOG_ERROR() << msg;								\
		handleStateMachineException(state, std::string("Error: ") + msg);	\
	}													\
	catch (...)											\
	{													\
		LOG_ERROR() << "Unknown error!";				\
		handleStateMachineException(state, std::string("Unknown error!"));	\
	}

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
const char BackEnd::ST_IDLE[] = "Idle";
const char BackEnd::ST_PREPARED[] = "Prepared";
const char BackEnd::ST_DESCRIBED[] = "Described";

const char BackEnd::ST_ACTIVE[] = "Active";
const char BackEnd::ST_READY[] = "Ready";
const char BackEnd::ST_RUNNING[] = "Running";
const char BackEnd::ST_PAUSED[] = "Paused";

const char BackEnd::EV_PREPARE[] = "Prepare";
const char BackEnd::EV_DESCRIBE[] = "Describe";
const char BackEnd::EV_UNDO[] = "Undo";
const char BackEnd::EV_CONFIGURE[] = "Configure";
const char BackEnd::EV_START[] = "Start";
const char BackEnd::EV_STOP[] = "Stop";
const char BackEnd::EV_PAUSE[] = "Pause";
const char BackEnd::EV_RESUME[] = "Resume";
const char BackEnd::EV_BREAK[] = "Break";
//_________________________________________________________________________________________________
/*
 * Constructs a BackEnd object, possibly launching one additional thread for an Ice server.
 */
BackEnd::BackEnd()
	: configManagerPtr(0), systemMgr_(*this), alarmLoggerAddress_(Default::alarmLoggerPortNum)
{
	// Initialize state machine status
	smErrorFlag = SM::NO_ERR;

	// Setup configuration manager
	setConfigManager();

	createEngine();
	engine.start();
}
//_________________________________________________________________________________________________
BackEnd::~BackEnd()
{
}
//_________________________________________________________________________________________________
void BackEnd::createAlarmLogger(const ::utl::net::SocketAddress & address)
{
	LOG_DEBUG() << "Creating alarm logger with address " << address;
	try
	{
		try
		{
			Server::find(address).addServant("AlarmLogger", new mt::AlarmLogger);
		}
		catch (const Server::NotCreatedException & e)
		{
			Server::create(address).addServant("AlarmLogger", new mt::AlarmLogger).start();
		}
		alarmLoggerAddress_ = address;
	}
	catch (const Ice::LocalException & e)
	{
		LOG_ERROR() << "Could not create alarm logger: " << e.ice_name();
	}
}
//_________________________________________________________________________________________________
mt::AlarmCallbackPrx BackEnd::alarmLogger()
{
	Ice::ObjectAdapterPtr adapter = ::mdaq::utl::Server::find(alarmLoggerAddress_).adapter();
	Ice::Identity id = adapter->getCommunicator()->stringToIdentity("AlarmLogger");

	Ice::ObjectPtr object = adapter->find(id);
	mt::AlarmCallbackPrx logger = 0;
	if (object)
	{
		logger = mt::AlarmCallbackPrx::checkedCast(adapter->createProxy(id));
	}
	return logger;
}
//_________________________________________________________________________________________________
const ConfigManager & BackEnd::configManager() const
{
	if (not configManagerPtr.get())
	{
		throw CmdException("No valid configuration manager has been setup!");
	}
	return *configManagerPtr.get();
}
//_________________________________________________________________________________________________
/** Sets up the manager of configurations available to the state machine.
 * @param repoUrl Path of the directory or URL of the database containing the configurations.
 */
void BackEnd::setConfigManager(const std::string & repoUrl)
{
	CCfg::Io::Url url(repoUrl);
	std::auto_ptr< ConfigManager> manager;
	if (url.isFile())
	{
		manager.reset(new ConfigFileManager(repoUrl));
	}
	else
	{
#if CFG_DB_SUPPORT_ENABLED
		manager.reset(new ConfigDatabaseManager(repoUrl));
#else
		LOG_ERROR() << "Could not create configuration manager for URL '" << repoUrl << "' (database support disabled).";
		return;
#endif
	}
	setConfigManager(manager);
}
//_________________________________________________________________________________________________
/** Sets up the manager of configurations available to the state machine.
 * @param manager New configuration manager.
 */
void BackEnd::setConfigManager(std::auto_ptr< ConfigManager > manager)
{
	configManagerPtr = manager;
}
//_________________________________________________________________________________________________
void BackEnd::loadPLLConfig(const float& ckwFreq_MHz)
{
	CoBoPLLManager(systemMgr_.getNode<CoBoNode>(currentHwNodeId_)).loadPLLConfig(ckwFreq_MHz);
}
//_________________________________________________________________________________________________
void BackEnd::writePLLRegister(const uint32_t& value)
{
	CoBoPLLManager(systemMgr_.getNode<CoBoNode>(currentHwNodeId_)).writePLLRegister(value);
}
//_________________________________________________________________________________________________
void BackEnd::resetPLL()
{
	CoBoPLLManager(systemMgr_.getNode<CoBoNode>(currentHwNodeId_)).resetPLL();
}
//_________________________________________________________________________________________________
/**
 * Creates the engine of the ECC state machine.
 */
void BackEnd::createEngine()
{

	dhsm::StateMachinePtr sm(new dhsm::StateMachine("SM_ECC"));
	sm->state(ST_IDLE);
	sm->state(ST_DESCRIBED);
	sm->state(ST_PREPARED);
	sm->state(ST_ACTIVE);
	sm->initial(ST_IDLE, dhsm::StateMachine::HISTORY);

	// Idle to Described
	sm->transition(EV_DESCRIBE, ST_IDLE, ST_DESCRIBED).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onDescribe);
	sm->transition(EV_UNDO, ST_DESCRIBED, ST_IDLE).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onUnDescribe);

	// Described to Prepared
	sm->transition(EV_PREPARE, ST_DESCRIBED, ST_PREPARED).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onPrepare);
	sm->transition(EV_UNDO, ST_PREPARED, ST_DESCRIBED).actions()
			.add(*this, &BackEnd::resetErrorFlag)
			.add(*this, &BackEnd::onUnPrepare);

	// Prepared to Active
	sm->transition(EV_CONFIGURE, ST_PREPARED, ST_ACTIVE).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onConfigure);
	sm->transition(EV_BREAK, ST_ACTIVE, ST_PREPARED).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onBreakup);

	// Active sub-machine
	dhsm::StateMachine & activeSm = sm->state(ST_ACTIVE).stateMachine("SM_ACTIVE");
	activeSm.initial(ST_READY, dhsm::StateMachine::NO_HISTORY);

	activeSm.transition(EV_START, ST_READY, ST_RUNNING).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onStart);
	activeSm.transition(EV_STOP, ST_RUNNING, ST_READY).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onStop);
	activeSm.transition(EV_STOP, ST_PAUSED, ST_READY).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onStop);
	activeSm.transition(EV_PAUSE, ST_RUNNING, ST_PAUSED).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onPause);
	activeSm.transition(EV_RESUME, ST_PAUSED, ST_RUNNING).actions()
		.add(*this, &BackEnd::resetErrorFlag)
		.add(*this, &BackEnd::onResume);

	engine.setStateMachine(sm);
}
//_________________________________________________________________________________________________
inline SM::State BackEnd::toSmState(const std::string & stateLabel)
{
	SM::State s = SM::OFF;
	if (stateLabel == ST_IDLE)
		s = SM::IDLE;
	else if (stateLabel == ST_DESCRIBED)
		s = SM::DESCRIBED;
	else if (stateLabel == ST_PREPARED)
		s = SM::PREPARED;
	else if (stateLabel == ST_ACTIVE or stateLabel == ST_READY)
		s = SM::READY;
	else if (stateLabel == ST_RUNNING)
		s = SM::RUNNING;
	else if (stateLabel == ST_PAUSED)
		s = SM::PAUSED;
	else
		LOG_ERROR() << "smState = " << stateLabel;

	return s;
}
//_________________________________________________________________________________________________
SM::State BackEnd::smState()
{
	return toSmState(engine.currentState().label());
}
//_________________________________________________________________________________________________
SM::Status BackEnd::smStatus()
{
	SM::Status s;
	s.error = smErrorFlag;
	s.state = smState();
	return s;
}
//__________________________________________________________________________________________________
void BackEnd::handleStateMachineException(SM::Error errorCode, const std::string & msg)
{
	smErrorFlag = errorCode;
	throw get::rc::SM::Exception(smStatus(), msg);
}
//__________________________________________________________________________________________________
/**
 * Returns status (state, transition, error) of ECC state machine.
 */
void BackEnd::getStateMachineStatus(SM::Status & response)
{
	response = smStatus();
}
//__________________________________________________________________________________________________
/**
 * Returns the list of data senders.
 */
 DataSenderSet BackEnd::dataSenders() const
{
	return systemMgr_.dataSenders();
}
//__________________________________________________________________________________________________
/**
 * Returns the mapping between CoBo or Mutant boards and their data routers.
 */
DataLinkSet BackEnd::dataLinks() const
{
	return systemMgr_.dataLinks();
}
//_________________________________________________________________________________________________
/**
 * Resets each read pointer to the beginning of its circular buffer.
 */
void BackEnd::resetCoBoReadPointers()
{
	connect("device", "ctrl");
	try
	{
		for (size_t asadIdx=0; asadIdx <4u; ++asadIdx)
		{
			std::string prefix = boost::lexical_cast< std::string >(asadIdx);
			writeRegister(std::string("readPtr_asad") + prefix, readRegister(std::string("memStart_asad") + prefix));
		}
	}
	catch (const CmdException &)
	{
		// Case of single AsAd firmware
		writeRegister(std::string("readPtr"), readRegister(std::string("memStart")));
	}
}
//_________________________________________________________________________________________________
/**
 * Starts run on current node.
 */
void BackEnd::daqStart()
{
	LOG_DEBUG() << ">>> daqStart";

	// Check run is not already started
	connect("device", "ctrl");
	if (readField("pipeCtrl", "intEnable") and readField("pipeCtrl", "firmwareEnable"))
	{
		LOG_WARN() << "Run has already been started!";
	}
	else
	{
		// Make sure firmware is not enabled
		writeField("pipeCtrl", "firmwareEnable", 0);

		// Make sure firmware does not generate interrupts
		writeField("pipeCtrl", "intEnable", 0);

		// Resets firmware fifo and state machine
		writeField("pipeCtrl", "reset", 1); // Resets write pointers to memStart values
		writeField("pipeCtrl", "backpressure", 1); // obsolete ?
		resetCoBoReadPointers(); // Resets read pointers to value of write pointers

		currentDaqCtrl()->daqStart();

		writeField("pipeCtrl", "backpressure", 0); // obsolete ?
		writeField("pipeCtrl", "reset", 0); // firmware is now in post-reset state

		// Firmware can now issue interrupts
		writeField("pipeCtrl", "intEnable", 1);
		// Go!
		writeField("pipeCtrl", "firmwareEnable", 1);
	}

	LOG_DEBUG() << "<<< daqStart";
}
//_________________________________________________________________________________________________
void BackEnd::daqStop()
{
	LOG_DEBUG() << ">>> daqStop";

	// Disable firmware acquisition
	connect("device", "ctrl");
	writeField("pipeCtrl", "firmwareEnable", 0);

	// Tell (interrupt) event consumer to stop after its current loop.
	currentDaqCtrl()->daqStop();
	writeField("pipeCtrl", "intEnable", 0);
	LOG_DEBUG() << "<<< daqStop";
}
//_________________________________________________________________________________________________
/**
 * Returns proxy to AsAd manager of current CoBo node.
 */
cobo::AsAdPulserMgrPrx BackEnd::asadPulserMgr()
{
	cobo::AsAdPulserMgrPrx mgr = 0;
	static const std::string errorMsg = "Could not connect to AsAd pulser manager. ";

	try
	{
		mgr = cobo::AsAdPulserMgrPrx::checkedCast(currentHwNodePrx());
	}
	catch (const Ice::LocalException & e)
	{
		throw CmdException(errorMsg + e.what());
	}
	if (not mgr)
	{
		throw CmdException(errorMsg);
	}
	return mgr;
}
//_________________________________________________________________________________________________
/**
 * Returns proxy to AsAd manager of current CoBo node.
 */
cobo::AsAdAlarmMonitorPrx BackEnd::asadAlarmMonitor()
{
	cobo::AsAdAlarmMonitorPrx mgr = 0;
	static const std::string errorMsg = "Could not connect to AsAd alarm monitor. ";

	try
	{
		mgr = cobo::AsAdAlarmMonitorPrx::checkedCast(currentHwNodePrx());
	}
	catch (const Ice::LocalException & e)
	{
		throw CmdException(errorMsg + e.what());
	}
	if (not mgr)
	{
		throw CmdException(errorMsg);
	}
	return mgr;
}
//_________________________________________________________________________________________________
/**
 * Returns proxy to alarm service of current CoBo node.
 */
mt::AlarmServicePrx BackEnd::alarmService()
{
	mt::AlarmServicePrx svce = 0;
	static const std::string errorMsg = "Could not connect to CoBo alarm service. ";

	try
	{
		svce = cobo::CtrlNodePrx::checkedCast(currentHwNodePrx())->getAlarmService();
	}
	catch (const Ice::LocalException & e)
	{
		throw CmdException(errorMsg + e.ice_name());
	}
	if (not svce)
	{
		throw CmdException(errorMsg);
	}
	return svce;
}
//__________________________________________________________________________________________________
/**
 * Subscribes ECC to alarms from current CoBo node.
 */
void BackEnd::subscribeAlarms(bool enabled)
{
	mt::AlarmServicePrx alarmSvce = alarmService();
	if (enabled)
	{
		alarmSvce->subscribe("ECC", alarmLogger());
	}
	else
	{
		alarmSvce->unsubscribe("ECC");
	}
}
//_________________________________________________________________________________________________
/**
 * Sets the ISPM multiplexer to debug the signals on AsAd.
 * See Fig. 14 and Tables 14-20 page 17 of "AsAd Firmware Description".
 * See "AsAd External Instruments Interface".
 * @param moduleName Name of device.
 * @param isp1SignalName Name of ISP1 signal.
 * @param isp2SignalName Name of ISP2 signal.
 */
void BackEnd::asadSetISPM(const std::string & moduleName, const std::string & isp1SignalName, const std::string & isp2SignalName)
{
	// Lookup register value.
	const uint64_t value = asad::utl::buildIspmRegisterValue(moduleName, isp1SignalName, isp2SignalName);
	LOG_DEBUG() << "Setting AsAd inspection manager register to " << std::showbase << std::hex << value << std::dec;

	// Set register
	connect("device", "asad");
	writeRegister("ISPM", value);
}
//_________________________________________________________________________________________________
void BackEnd::connect(const std::string & objectType, const std::string & objectName, bool isBatched)
{
	std::string daqFlowType = objectType;
	ba::to_lower(daqFlowType);
	if (ba::istarts_with(daqFlowType, std::string("daq-")))
	{
		daqFlowType = daqFlowType.substr(4);
	}

	if (daqFlowType == "ice")
	{
		connectRunProcessor(objectName, "ICE");
	}
	else if (daqFlowType == "tcp")
	{
		connectRunProcessor(objectName, "TCP");
	}
	else if (daqFlowType == "zbuf")
	{
		connectRunProcessor(objectName, "ZBUF");
	}
	else if (daqFlowType == "narval" or daqFlowType == "fdt")
	{
		connectRunProcessor(objectName, "FDT");
	}
	else if (daqFlowType == "debug")
	{
		connectRunProcessor(objectName, "Debug");
	}
	else if (daqFlowType == "debug-frames")
	{
		connectRunProcessor(objectName, "DebugFrames");
	}
	else if (daqFlowType == "none")
	{
		disconnectRunProcessor();
	}
	else
	{
		EccBackEnd::connect(objectType, objectName, isBatched);
	}
}
//_________________________________________________________________________________________________
/**
 * Returns all configuration IDs that are available to the ECC.
 * @param[out] configIds A set of configuration IDs.
 */
void BackEnd::getConfigIds(ConfigIdSet & configIds)
{
	try
	{
		configManager().getConfigIds(configIds);
	}
	catch (const std::exception & e)
	{
		handleStateMachineException(SM::NO_ERR, e.what());
	}
}
//_________________________________________________________________________________________________
/**
 * Returns all configuration IDs of given type.
 * @param subConfigType Type of configuration.
 * @param[out] subConfigIds A set of configuration identifiers.
 */
void BackEnd::getSubConfigIds(const std::string & subConfigType, SubConfigIdSet & subConfigIds) const
{
	try
	{
		configManager().getSubConfigIds(ConfigId::stringToSubType(subConfigType), subConfigIds);
	}
	catch (const std::exception & e)
	{
		throw CmdException(e.what());
	}
}
//__________________________________________________________________________________________________
std::string BackEnd::addSubConfig(const std::string & subConfigType, const std::string & subConfig, const std::string & subConfigId) const
{
	return configManager().addSubConfig(ConfigId::stringToSubType(subConfigType), subConfig, subConfigId);
}
//__________________________________________________________________________________________________
std::string BackEnd::getSubConfig(const std::string & subConfigType, const std::string & subConfigId) const
{
	std::string xmlConfig;
	configManager().getSubConfig(ConfigId::stringToSubType(subConfigType), subConfigId, xmlConfig);
	return xmlConfig;
}
//__________________________________________________________________________________________________
std::string BackEnd::getConfigSet(const ConfigId & configId) const
{
	return configManager().getConfigSet(configId);
}
//__________________________________________________________________________________________________
/**
 * Returns the contents of the last configuration used as XML.
 */
std::string BackEnd::getLastConfigSet() const
{
	return getConfigSet(lastConfigId);
}
//__________________________________________________________________________________________________
/**
 * Describes the current hardware node, based on a configuration ID or a full configuration URL.
 */
void BackEnd::loadHwDescription(const std::string & descrRef)
{
	if (boost::starts_with(descrRef, "<?xml"))
	{
		std::string xmlString = descrRef;
		try
		{
			CCfg::Io::BufferStorage buffer;
			buffer.write(xmlString);
			CCfg::Io::XmlFormat xml;
			std::auto_ptr<Object> descrObj = xml.parse(buffer, CCfg::Io::NoOpIncludeHandler::getInstance());
			mdaq::EccBackEnd::loadHwDescription(CCfg::CConfig(descrObj.get()));
		}
		catch (const std::exception & e)
		{
			throw mdaq::utl::CmdException(e.what());
		}
	}
	else
	{
		// Try hardware config ID
		try
		{
			std::auto_ptr< CCfg::Object > descrObj = configManager().getHardwareDescription(descrRef);
			mdaq::EccBackEnd::loadHwDescription(CCfg::CConfig(descrObj.get()));
		}
		catch (const ConfigManager::ConfigIdNotFound &)
		{
			try
			{
				// Try configuration file
				CCfg::Io::Url url(descrRef);
				if (url.isFile())
				{
					mdaq::EccBackEnd::loadHwDescription(descrRef);
				}
				// Try configuration URL
				else
				{
					std::auto_ptr< CCfg::Object > descrObj = configManager().getConfig(url);
					mdaq::EccBackEnd::loadHwDescription(CCfg::CConfig(descrObj.get()));
				}
			}
			catch (const ConfigManager::ConfigIdNotFound &)
			{
				throw CmdException(std::string("Could not find hardware description '") + descrRef + '\'');
			}
			catch (const CCfg::Exception & e)
			{
				throw CmdException(std::string("Could not load hardware description '") + descrRef + '\'');
			}
			catch (const CmdException & e)
			{
				throw CmdException(std::string("Could not load hardware description '") + descrRef + "\': " + e.reason);
			}
		}
	}
}
//__________________________________________________________________________________________________
/**
 * Sets up hardware topology (which boards are in the system...).
 */
void BackEnd::describe(const std::string & subConfigId)
{
	dhsm::Message<std::string> msg(EV_DESCRIBE, subConfigId);
	engine.step(msg);
}
//_________________________________________________________________________________________________
void BackEnd::removeAllHwNodes()
{
	systemMgr_.removeNodes();
	EccBackEnd::removeAllHwNodes();
}
//__________________________________________________________________________________________________
/**
 * Describes current node.
 */
void BackEnd::describeNode(const std::string nodeType, const std::string & instanceId, const std::string & cfgUrl)
{
	// Load description
	loadHwDescription(cfgUrl);

	// Rename node
	EccBackEnd::renameHwNode(buildNodeId(nodeType, instanceId));

	// Register CoBo boards or MuTanT module
	systemMgr_.addNewNode(nodeType, instanceId, currentHwNode());
}
//__________________________________________________________________________________________________
/**
 * Sets up hardware topology (which boards are in the system...).
 */
void BackEnd::onDescribe(std::string subConfigId)
{
	LOG_INFO() << "Describe: loading topology '" << subConfigId << "'";
	try
	{
		// Store config ID
		lastConfigId.clear();
		lastConfigId.erase(ConfigId::Describe);
		lastConfigId.insert(ConfigId::Describe, subConfigId);

		// Get description configuration
		std::auto_ptr< CCfg::Object > descriptionObj = configManager().getDescribeConfig(subConfigId);
		CCfg::CConfig config = descriptionObj.get();

		// Remove all nodes
		LOG_DEBUG() << "Deleting already existing nodes...";
		systemMgr_.removeNodes();
		removeAllHwNodes();

		// Loop over node types (CoBo, Mutant, ...)
		CCfg::View::InheritersIterator inheritersIter =
				config("Node", "*")("Instance", "*").getInheriterAttributes();
		while (inheritersIter.hasNext())
		{
			// Do not forget to call next before accessing the first element.
			CCfg::Object & nodeObj = inheritersIter.next()->toObject();
			// Skip node templates
			if (nodeObj.isDefaultConfig())
				continue;

			// Check node instance is active
			if (not isNodeInstanceActive(nodeObj))
				continue;

			// Get type of node ('CoBo', 'Mutant', ...)
			const std::string nodeType = nodeObj.getParent()->getIndex();
			// Get name of node
			const std::string instanceId = nodeObj.getIndex();
			// Get address
			::mdaq::utl::net::SocketAddress endpoint;
			endpoint.fromCConfig(config("Node", nodeType)("Instance", instanceId)("Endpoint"));

			// Check whether batched processing should be enabled
			bool isBatched = config("Node", nodeType)("Instance", instanceId)("batchProcessing").getBoolValue(false);
			if (nodeType == MutantNode::NODE_TYPE)
			{
				isBatched = false;
			}

			// Connect to node
			LOG_DEBUG() << "Connecting to '" << nodeType << "' node '" << instanceId << "' at " << endpoint;
			try
			{
				connect("node", endpoint.toString(), isBatched);
			}
			catch (const CmdException & e)
			{
				throw CmdException(std::string("Could not connect to '") + nodeType
						+ "' node '" + instanceId + "' at " + endpoint.toString() + ". " + e.reason);
			}

			// Get ID of hardware description
			const std::string hardwareId = config("Node", nodeType)("Instance", instanceId)("HardwareConfigId").find().getStringValue();

			// Describe node
			describeNode(nodeType, instanceId, hardwareId);
		}

		// Set the mutantConfig register of the CoBo boards.
		systemMgr_.describe();

		//std::cout << "Data senders: " << dataSenders() << std::endl;

	}
	CATCH_SM_EXCEPTIONS(SM::WHEN_DESCRIBE)
}
//__________________________________________________________________________________________________
std::string BackEnd::buildNodeId(const std::string & nodeType, const std::string & instanceId)
{
	return nodeType + "[" + instanceId + "]";
}
//__________________________________________________________________________________________________
/**
 * Checks node instance is marked as active in description file.
 */
bool BackEnd::isNodeInstanceActive(CCfg::Object & instance)
{
	bool isActive = true;
	try
	{
		isActive = CCfg::CConfig(&instance)("isActive").getBoolValue();
	}
	catch (const CCfg::Exception & e)
	{
		LOG_WARN() << e.what();
	}
	return isActive;
}
//__________________________________________________________________________________________________
/**
 * Inter-board communications, clock synchronizations, calibrations
 */
void BackEnd::prepare(const std::string & subConfigId)
{
	dhsm::Message<std::string> msg(EV_PREPARE, subConfigId);
	engine.step(msg);
}
//__________________________________________________________________________________________________
/**
 * Prepares a node.
 * @param nodeType Type of the node.
 * @param instanceId Which instance to prepare among nodes of this type.
 * @param configRef URL of ID of the configuration.
 */
void BackEnd::prepareNode(const std::string & nodeType, const std::string & instanceId, const std::string & configRef)
{
	// Get configuration
	std::auto_ptr< CCfg::Object > rootObj;
	// Try configuration ID
	try
	{
		rootObj = configManager().getPrepareConfig(configRef);
	}
	catch (const ConfigManager::ConfigIdNotFound &)
	{
		// Try configuration file or URL
		try
		{
			CCfg::Io::Url url(configRef);
			rootObj = configManager().getConfig(url);
		}
		catch (const ConfigManager::ConfigIdNotFound &)
		{
			throw CmdException(std::string("Could not find configuration '") + configRef + "' to prepare node.");
		}
	}

	// Handle structure of old GetController configuration files
	CCfg::CConfig config = rootObj.get();
	if (not CCfg::CConfig(rootObj.get())("isActive").exists()
			and CCfg::CConfig(rootObj.get())("CoBo", instanceId)("isActive").exists())
	{
		config = CCfg::CConfig(rootObj.get())("CoBo", instanceId);
	}

	// Apply
	prepareNode(nodeType, instanceId, config);
}
//__________________________________________________________________________________________________
/**
 * Prepares a node.
 * @param nodeType Type of the node.
 * @param instanceId Which instance to prepare among nodes of this type.
 * @param config Configuration.
 */
void BackEnd::prepareNode(const std::string & nodeType, const std::string & instanceId, const CCfg::CConfig & config)
{
	systemMgr_.getNode(nodeType, instanceId).call_prepare(config);
}
//__________________________________________________________________________________________________
/**
 * Inter-board communications, clock synchronizations, calibrations.
 * @param subConfigId Identifier of the configuration to use.
 * Prepares all nodes:
 *   - prepare each CoBo board
 *   - prepares each Mutant module
 *   - find out which slots the CoBo boards are in, using Mutant's alignment procedure
 */
void BackEnd::onPrepare(std::string subConfigId)
{
	LOG_INFO() << "Prepare: using configuration '" << subConfigId << "' for inter-board communications, clock synchronizations, calibrations...";
	try
	{
		// Store config ID
		lastConfigId.erase(ConfigId::Prepare);
		lastConfigId.insert(ConfigId::Prepare, subConfigId);

		// Get configuration
		std::auto_ptr< CCfg::Object > cfgObj = configManager().getPrepareConfig(subConfigId);
		CCfg::CConfig config = cfgObj.get();

		// Prepare system
		systemMgr_.prepare(config);
	}
	CATCH_SM_EXCEPTIONS(SM::WHEN_PREPARE)
}
//__________________________________________________________________________________________________
/**
 * Triggers the "Undo" transition.
 */
void BackEnd::reset()
{
	engine.step(EV_UNDO);
}
//__________________________________________________________________________________________________
/**
 * Forgets system topology.
 */
void BackEnd::onUnDescribe()
{
	try
	{
		// Remove all nodes
		systemMgr_.removeNodes();
		removeAllHwNodes();
	}
	CATCH_SM_EXCEPTIONS(SM::WHEN_RESET)
}
//__________________________________________________________________________________________________
/**
 * Forgets system topology.
 */
void BackEnd::onUnPrepare()
{
	LOG_INFO() << "Un-preparing system";
	systemMgr_.unprepare();
}
//__________________________________________________________________________________________________
/**
 * Triggers the "Configure" transition.
 */
void BackEnd::configure(const std::string & subConfigId, const std::string & dataRouters)
{
	if (engine.currentState().label() == ST_PREPARED)
	{
		dhsm::Message< ConfigParams > msg(EV_CONFIGURE, std::make_pair(subConfigId, dataRouters));
		engine.step(msg);
	}
}
//__________________________________________________________________________________________________
void BackEnd::configureNode(const std::string & nodeType, const std::string & instanceId, const std::string & configRef)
{
	// Get configuration
	std::auto_ptr< CCfg::Object > cfgObj;
	// Try configuration ID
	try
	{
		cfgObj = configManager().getConfigureConfig(configRef);
	}
	catch (const ConfigManager::ConfigIdNotFound &)
	{
		// Try configuration file or URL
		try
		{
			CCfg::Io::Url url(configRef);
			cfgObj = configManager().getConfig(url);
		}
		catch (const ConfigManager::ConfigIdNotFound &)
		{
			throw CmdException(std::string("Could not find configuration '") + configRef + "' to configure node.");
		}
	}

	// Apply
	CCfg::CConfig config = cfgObj.get();
	configureNode(nodeType, instanceId, config);
}
//__________________________________________________________________________________________________
/**
 * Configures a node.
 * @param nodeType Type of the node.
 * @param instanceId Which instance to configure among nodes of this type.
 * @param config Configuration.
 */
void BackEnd::configureNode(const std::string & nodeType, const std::string & instanceId, const CCfg::CConfig & nodeCfg)
{
	systemMgr_.getNode(nodeType, instanceId).call_configure(nodeCfg);
}
//__________________________________________________________________________________________________
/**
 * Loads physics parameters into the registers of the boards (thresholds, ...).
 */
void BackEnd::onConfigure(ConfigParams parameters)
{
	try
	{
		const std::string subConfigId = parameters.first;
		LOG_INFO() << "Configure: loading physics parameters from configuration '" << subConfigId << "'";

		// Store config ID
		lastConfigId.erase(ConfigId::Configure);
		lastConfigId.insert(ConfigId::Configure, subConfigId);

		// Get configuration
		std::auto_ptr<CCfg::Object> descriptionObj = configManager().getConfigureConfig(subConfigId);
		CCfg::CConfig config = descriptionObj.get();

		// Configure system
		systemMgr_.configure(config, parameters.second);

	}
	CATCH_SM_EXCEPTIONS(SM::WHEN_CONFIGURE)
}
//__________________________________________________________________________________________________
/**
 * Triggers the 'Start' transition.
 */
void BackEnd::start()
{
	engine.step(EV_START);
}
//__________________________________________________________________________________________________
/**
 * Starts the run, i.e. allows hardware to deliver data.
 */
void BackEnd::onStart()
{
	try
	{
		try
		{
			systemMgr_.start();
		}
		catch (const Ice::LocalException & e)
		{
			throw CmdException(std::string("Could not start run (") + e.ice_name() + ')');
		}
	}
	CATCH_SM_EXCEPTIONS(SM::WHEN_START)
}
//__________________________________________________________________________________________________
/**
 * Triggers the 'Stop' transition.
 */
void BackEnd::stop()
{
	engine.step(EV_STOP);
}
//__________________________________________________________________________________________________
/**
 * Ends the run, preventing hardware from delivering data.
 * Leads to Ready state.
 */
void BackEnd::onStop()
{
	try
	{
		try
		{
			systemMgr_.stop();
		}
		catch (const Ice::LocalException & e)
		{
			throw CmdException(std::string("Could not stop run (") + e.ice_name() + ')');
		}
	}
	CATCH_SM_EXCEPTIONS(SM::WHEN_STOP)
}
//__________________________________________________________________________________________________
/**
 * Triggers the 'Pause' transition.
 */
void BackEnd::pause()
{
	engine.step(EV_PAUSE);
}
//__________________________________________________________________________________________________
/**
 * Pauses the run, i.e. forbids hardware to deliver data.
 * Leads to Paused state.
 */
void BackEnd::onPause  ()
{
	try
	{
		systemMgr_.pause();
	}
	CATCH_SM_EXCEPTIONS(SM::WHEN_PAUSE)
}
//__________________________________________________________________________________________________
/**
 * Triggers the 'Resume' transition.
 */
void BackEnd::resume()
{
	engine.step(EV_RESUME);
}
//__________________________________________________________________________________________________
/**
 * Resumes the run, i.e. allows again hardware to deliver data.
 */
void BackEnd::onResume ()
{
	try
	{
		systemMgr_.resume();
	}
	CATCH_SM_EXCEPTIONS(SM::WHEN_RESUME)
}
//__________________________________________________________________________________________________
/**
 * Triggers the 'Breakup' transition.
 */
void BackEnd::breakup()
{
	engine.step(EV_BREAK);
}
//__________________________________________________________________________________________________
/**
 * Invalidates physics parameters setup in the boards.
 * Leads to Idle state (so physics parameters will have to be setup again).
 */
void BackEnd::onBreakup()
{
	try
	{
		onStop();
		systemMgr_.breakup();
	}
	CATCH_SM_EXCEPTIONS(SM::WHEN_BREAKUP)
}
//__________________________________________________________________________________________________
void BackEnd::resetErrorFlag()
{
	smErrorFlag = SM::NO_ERR;
}
//__________________________________________________________________________________________________
std::string BackEnd::specialCmd(const std::string & cmdLine)
{
	std::string ret = "";

	std::istringstream cmdArgs(cmdLine);
	std::string cmdName;
	cmdArgs >> cmdName;
	LOG_INFO() << "Running special command '" << cmdName << "'";
	boost::algorithm::to_lower(cmdName);

	try
	{
		if ("asad-inject" == cmdName)
		{
			long amplitude = 0;
			cmdArgs >> amplitude;

			asadPulserMgr()->triggerPulser(amplitude);
		}
		else if ("connect" == cmdName)
		{
			std::string objectType;
			std::string objectName;

			cmdArgs >> objectType;
			cmdArgs >> objectName;

			connect(objectType, objectName);
		}
		else if ("write" == cmdName)
		{
			std::string registerName, fieldName;
			uint64_t value;
			cmdArgs >> registerName;
			cmdArgs >> fieldName;
			cmdArgs.unsetf(std::ios::dec);
			cmdArgs >> value;

			writeField(registerName, fieldName, value);
		}
		else if ("writereg" == cmdName)
		{
			std::string registerName;
			uint64_t value;
			cmdArgs >> registerName;
			cmdArgs.unsetf(std::ios::dec);
			cmdArgs >> value;
			if (cmdArgs.fail())
			{
				throw get::rc::SM::Exception(smStatus(), std::string("Could not parse arguments for command '") + cmdName + "'.");
			}
		}
		else
		{
			throw get::rc::SM::Exception(smStatus(), std::string("Special command '") + cmdName + "' is unknown.");
		}
	}
	catch (const CmdException & e)
	{
		throw get::rc::SM::Exception(smStatus(), e.reason);
	}

	return ret;
}
//__________________________________________________________________________________________________
/**
 * Sets writing clock CKW frequency of current CoBo node.
 *   - sets the LCM1 register of currently selected AsAd boards for given frequency
 *   - sets the CoBo 'ckwDivider' register for the CKW frequency divider: CKW = 100 MHz / 2**ckdDivider
 *   - configures M-CoBo PLL device for given frequency.
 *
 * @param frequency_MHz Frequency [MHz].
 * @param pllConfigEnabled Debugging flag to disable configuration of the CoBo PLL device.
 * @param regs Map of non-default PLL register values.
 */
void BackEnd::coboSetWritingClockFrequency(const float & frequency_MHz, bool pllConfigEnabled, const cobo::PllRegisterMap & regs)
{
	systemMgr_.getCoBoNode(currentHwNodeId_).setCoBoWritingClockFrequency(frequency_MHz, pllConfigEnabled, regs);
}
//__________________________________________________________________________________________________
/**
 * Sets reading clock CKR frequency of current CoBo node.
 *  - sets the CoBo 'ckrDivider' register for the CKR frequency divider: CKR = 100 MHz / 2**ckdDivider
 * @param frequency_MHz Frequency [MHz].
 * The only frequency value allowed is 25 MHz.
 */
void BackEnd::coboSetReadingClockFrequency(const float & frequency_MHz)
{
	systemMgr_.getCoBoNode(currentHwNodeId_).setCoBoReadingClockFrequency(frequency_MHz);
}
//__________________________________________________________________________________________________
/**
 * Switches current CoBo node to slow or fast control mode.
 * @param slowMode Whether to switch from fast to slow (true) or from slow to fast control mode.
 */
void BackEnd::coboSwitchToSlowControlMode(bool slowMode)
{
	// Store mask of active AGET chips
	// const uint16_t agetMask = ctrlDevice->readField("aget_sc_configuration", "SC_En");
	CoBoNode& coboNode = systemMgr_.getNode<CoBoNode>(currentHwNodeId_);

	LOG_DEBUG() << "Switching to CoBo " << (slowMode ? "slow" : "fast") << " control mode";
	if (slowMode)
	{
		coboNode.switchToSlowControl();
	}
	else
	{
		// Switch to hit register mode
		coboNode.switchToFastControl();
		// Restart firmware state machine
		coboNode.hwNode().regProc("ctrl").write("pipeCtrl", "firmwareEnable", 1);
	}
}
//__________________________________________________________________________________________________
/**
 * Selects single channel (per chip) of current CoBo node to always read.
 */
void BackEnd::coboSelectSingleReadoutChannel(const size_t& chanId)
{
	// Store mask of active AGET chips
	// const uint16_t agetMask = ctrlDevice->readField("aget_sc_configuration", "SC_En");

	LOG_DEBUG() << "Selecting single channel no. " << chanId << " to read";
	systemMgr_.getCoBoNode(currentHwNodeId_).selectSingleAgetReadoutChannel(chanId);
}
//__________________________________________________________________________________________________
/**
 * Selects single channel (per chip) of current CoBo node for AGET test mode.
 */
void BackEnd::coboSelectSingleTestChannel(const size_t & chanId)
{
	LOG_DEBUG() << "Selecting single channel no. " << chanId << " to test";
	systemMgr_.getCoBoNode(currentHwNodeId_).selectSingleAgetTestChannel(chanId);
}
//__________________________________________________________________________________________________
/**
 * Modifies LSB threshold of a single channel (per chip) of current CoBo node for test mode.
 * @param chanId Index of the single test channel to read (0 to 67). If 68, set all.
 * @param threshold Least significant bits of channel threshold (from 0 to 15).
 */
void BackEnd::coboSetAgetChannelLsbThreshold(const size_t & chanId, const uint8_t & threshold)
{
	systemMgr_.getCoBoNode(currentHwNodeId_).setAgetLsbThreshold(chanId, threshold);
}
//__________________________________________________________________________________________________
/**
 * 	Sets up CoBo-Mutant link mode:
 *	  - Normal Operation (Input buffer enabled)
 *	  - CoBo Disconnected from Mutant network
 *	  - GMC Feedback Mode (Alignment)
 *	  - WSCA Feedback Mode (Alignment)
 */
void BackEnd::setupCoBoMutantLink(get::cobo::MutantLinkMode mode)
{
	try
	{
		const std::string proxyString = currentHwNode().endpoint.iceProxyString("MutantLinkManager", "default", proxyOptions);
		Ice::ObjectPrx base = currentHwNodePrx()->ice_getCommunicator()->stringToProxy(proxyString);
		cobo::MutantLinkManagerPrx linkManager = cobo::MutantLinkManagerPrx::checkedCast(base);
		linkManager->setupMutantLink(mode);
	}
	catch (const Ice::LocalException & e)
	{
		throw CmdException(std::string("Could not connect to CoBo Mutant link manager: ") + e.ice_name());
	}
}
//__________________________________________________________________________________________________
/**
 * Get proxy for alignment server of current Mutant node.
 */
mutant::AlignmentServerPrx BackEnd::mutantAlignmentServer()
{
	// Get proxy  for alignment server
	const std::string proxyString = currentHwNode().endpoint.iceProxyString("MutantAlignmentServer", "default", proxyOptions);
	Ice::ObjectPrx base = currentHwNodePrx()->ice_getCommunicator()->stringToProxy(proxyString);
	mutant::AlignmentServerPrx alignServer = mutant::AlignmentServerPrx::checkedCast(base);
	return alignServer;
}
//_________________________________________________________________________________________________
/**
 * Returns proxy to LED manager of current CoBo node.
 */
cobo::LedManagerPrx BackEnd::ledManager()
{
	cobo::LedManagerPrx svce = 0;
	static const std::string errorMsg = "Could not connect to CoBo LED manager. ";

	try
	{
		svce = cobo::LedManagerPrx::checkedCast(currentHwNodePrx());
	}
	catch (const Ice::LocalException & e)
	{
		throw CmdException(errorMsg + e.ice_name());
	}
	if (not svce)
	{
		throw CmdException(errorMsg);
	}
	return svce;
}
//__________________________________________________________________________________________________
} // namespace rc
} // namespace get
