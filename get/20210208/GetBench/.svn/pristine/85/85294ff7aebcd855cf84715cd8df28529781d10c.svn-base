/**
 * @file CoBoNode.cpp
 * @date 15 oct. 2012
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

#include "get/rc/CoBoNode.h"
#include "get/asad/utl.h"
#include "get/rc/BackEnd.h"
#include "get/rc/DataLinkSet.h"
#include "get/rc/CoBoPLLManager.h"
#include "mdaq/rc/HardwareNode.h"
#include "mdaq/hw/Control.h"
#include "mdaq/utl/CmdException.h"
#include "utl/Logging.h"
#include <IceUtil/Thread.h>
#include ICE_HEADER(LocalException.h)
#include <boost/cstdint.hpp>
#ifndef UINT64_C
#define UINT64_C(literal) literal ## ull
#endif
#include <boost/lexical_cast.hpp>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <string>

using ::utl::BitFieldHelper;
using ::mdaq::utl::CmdException;
using namespace mdaq::hw;

namespace get
{
namespace rc
{

const uint16_t CoBoNode::KNOWN_AGET_VERSIONS[NUM_KNOWN_AGET_VERSIONS] = {0x201, 0x202, 0x203}; ///< ASTRE = 0x203
const std::string CoBoNode::NODE_TYPE("CoBo");

// --------------------------------------------------------------------------------------------------------------------
CoBoNode::CoBoNode(const std::string& coboInstId, mdaq::rc::HardwareNode* _hwNode, get::rc::BackEnd* backend)
	: Node(NodeId(NODE_TYPE, coboInstId), backend, _hwNode), agetMask_(0), asadMask_(0), isActive_(true), useExternalTimeReset_(true)
{
	if (_hwNode)
	{
		NODE_DEBUG() << "Creating manager for board \'" << id().toString() << '\'';
	}
}
// --------------------------------------------------------------------------------------------------------------------
CoBoNode::~CoBoNode()
{
	if (hwNode_)
	{
		// Set the 'initDone' bit of CoBo's 'mutantConfig' Mutant configuration register to 0.
		setInitDone(false);
	}
}
// --------------------------------------------------------------------------------------------------------------------
DataSenderId CoBoNode::dataSenderId() const
{
	return DataSenderId(id().toString());
}
// --------------------------------------------------------------------------------------------------------------------
DataSenderSet CoBoNode::dataSenders() const
{
	DataSenderSet senders;
	senders.insert(dataSenderId());
	return senders;
}
// --------------------------------------------------------------------------------------------------------------------
get::cobo::CtrlNodePrx CoBoNode::ctrlPrx()
{
	return get::cobo::CtrlNodePrx::checkedCast(hwNode().hwServerCtrlPrx);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Returns proxy for AsAd pulser manager in CoBo.
 */
get::cobo::AsAdPulserMgrPrx CoBoNode::asadPulserMgr()
{
	return get::cobo::AsAdPulserMgrPrx::checkedCast(hwNode().hwServerCtrlPrx);
}

// --------------------------------------------------------------------------------------------------------------------
/**
 * Returns proxy for CoBo monitor of AsAd alarms.
 */
get::cobo::AsAdAlarmMonitorPrx CoBoNode::asadAlarmMonitor()
{
	return get::cobo::AsAdAlarmMonitorPrx::checkedCast(hwNode().hwServerCtrlPrx);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Returns proxy for CoBo LED manager.
 */
cobo::LedManagerPrx CoBoNode::ledManager()
{
	return get::cobo::LedManagerPrx::checkedCast(hwNode().hwServerCtrlPrx);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Returns proxy for CoBo manager of Mutant link.
 */
cobo::MutantLinkManagerPrx CoBoNode::mutantLinkManager()
{
	std::string proxyOptions;
#if ICE_INT_VERSION >= 30500
	proxyOptions = "-e 1.0";
#endif
	const std::string proxyString = hwNode().endpoint.iceProxyString("MutantLinkManager", "default", proxyOptions);
	Ice::ObjectPrx base = hwNode().hwServerCtrlPrx->ice_getCommunicator()->stringToProxy(proxyString);
	return cobo::MutantLinkManagerPrx::checkedCast(base);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Returns proxy for 'ctrl' device.
 */
mdaq::hw::DevicePrx& CoBoNode::ctrl()
{
	return hwNode().findDevice("ctrl");
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Returns proxy for 'asad' slow control device.
 */
mdaq::hw::DevicePrx& CoBoNode::asad()
{
	return hwNode().findDevice("asad");
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Returns proxy for 'aget' slow control device.
 */
mdaq::hw::DevicePrx& CoBoNode::aget()
{
	return hwNode().findDevice("aget");
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Returns proxy for zero suppression device for a given AsAd board.
 * @param asadIdx AsAd board ID
 */
mdaq::hw::DevicePrx& CoBoNode::zeroSuppressDevice(size_t asadIdx)
{
	std::ostringstream deviceName;
	deviceName << "zeroSuppress";
	if (asadIdx > 0)
	{
		deviceName << "_AsAd" << asadIdx;
	}
	return hwNode().findDevice(deviceName.str());
}
// --------------------------------------------------------------------------------------------------------------------
void CoBoNode::sleep_ms(size_t delay_ms)
{
	IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(delay_ms));
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Creates a mask of 16 bits indicating with of the 16 AGET ASICs are active.
 */
void CoBoNode::setAgetMask(CCfg::CConfig nodeCfg)
{
	isActive_ = nodeCfg("isActive").getBoolValue();
	asadMask_ = 0;
	agetMask_ = 0;

	if (isActive_)
	{
		for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
		{
			CCfg::CConfig asadCfg = nodeCfg("AsAd", asadIdx);
			bool isBoardActive = asadCfg("Control")("isActive").getBoolValue();
			if (isBoardActive)
			{
				::utl::BitFieldHelper<uint16_t>::setField(asadMask_, asadIdx, 1, 1);
				// Loop over 4 AGET chips
				for (size_t agetIdx=0; agetIdx < 4; ++agetIdx)
				{
					bool isChipActive = asadCfg("Aget", agetIdx)("Control")("isActive").getBoolValue();
					::utl::BitFieldHelper<uint16_t>::setField(agetMask_, 4*asadIdx+agetIdx, 1, (isChipActive?1u:0));
				}
			}
		}
	}
	NODE_DEBUG() << "AsAd mask: " << std::showbase << std::hex << asadMask_ << std::dec
				<< " ; AGET mask: " << std::showbase << std::hex << agetMask_ << std::dec;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Tests whether ASIC is active.
 */
bool CoBoNode::isAgetActive(size_t asadId, size_t agetId) const
{
	return utl::BitFieldHelper< uint16_t >::getField(agetMask_, 4*asadId + agetId, 1);
}

/**
 * Read AGET version register.
 */
uint16_t CoBoNode::agetVersion(size_t asadId, size_t agetId)
{
	// Select AGET to read from
	ctrl()->writeField("agetChipSelect", "d", 4 * asadId + agetId);
	// Read AGET register 5
	return aget()->readRegister("reg5");
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Returns ID (0 to 15) of first active chip
 */
size_t CoBoNode::firstActiveChip() const
{
	size_t id=0;
	for (size_t id=0; id<16; ++id)
	{
		if (BitFieldHelper< uint16_t >::getBit(agetMask_, id))
			break;
	}
	return id;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Tests whether AsAd is active.
 */
bool CoBoNode::isAsadActive(size_t asadId) const
{
	return utl::BitFieldHelper< uint16_t >::getField(asadMask_, asadId, 1) != 0;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Tests whether CoBo is active.
 */
bool CoBoNode::isCoBoActive() const
{
	return isActive_;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Checks whether AsAd board is plugged in into the CoBo board.
 * Reads CoBo read-only register giving status of the AsAd boards.
 * @param asadIdx Index of the AsAd board.
 */
bool CoBoNode::isAsAdConnected(size_t asadIdx)
{
	const std::string fieldName = std::string("PLG") + boost::lexical_cast< std::string >(asadIdx);
	bool isConnected = (1 == ctrl()->readField("asadConnection", fieldName));
	return isConnected;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Read AsAd board AL bit.
 * Its meaning varies depending on PLG bit.
 * @param asadIdx Index of the AsAd board.
 * @see GET-AS-002-0004
 */
uint8_t CoBoNode::asadAlarmBit(const size_t & asadIdx)
{
	return isAsAdPowerSupplied(asadIdx);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Checks whether AsAd board is connected to its power supply.
 * Reads AsAd alarm bits.
 * @param asadIdx Index of the AsAd board.
 * @see GET-AS-002-0004
 */
bool CoBoNode::isAsAdPowerSupplied(size_t asadIdx)
{
	const uint8_t alarmBits = ctrl()->readField("asadStatus", "alarm");
	const bool isPowerSupplied = utl::BitFieldHelper< uint8_t >::getBit(alarmBits, asadIdx);
	return isPowerSupplied;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Checks whether AsAd board was powered on.
 * @see GET-AS-002-0004
 */
bool CoBoNode::isAsAdPoweredOn(const size_t & asadIdx)
{
	const std::string fieldName = std::string("powerON") + boost::lexical_cast< std::string >(asadIdx);
	bool pwu = (1 == ctrl()->readField("asadEnable", fieldName));
	return pwu;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Selects AsAd board(s) for slow control.
 */
void CoBoNode::setAsAdSlowControlMask(uint16_t mask)
{
	try
	{
		ctrl()->writeField("ASAD_sc_mask", "select", mask);
	}
	catch (const CmdException & e)
	{
		NODE_WARN() << e.reason;
		NODE_WARN() << "You need to update your hardware description file!";
		throw;
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Selects AsAd board(s) to read from.
 */
void CoBoNode::selectAsAdForRead(size_t asadId)
{
	ctrl()->writeField("asadSelect", "d", asadId);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Selects AGET chips for slow control.
 */
void CoBoNode::setAgetSlowControlMask(uint16_t mask)
{
	ctrl()->writeField("aget_sc_configuration", "SC_En", mask);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Selects a unique AsAd board for a slow control command.
 * @param asadIdx Index of the AsAd board to select.
 */
void CoBoNode::selectAsAdForSlowControl(size_t asadIdx)
{
	const uint16_t mask = (1ul << asadIdx);
	setAsAdSlowControlMask(mask);
	selectAsAdForRead(asadIdx);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Selects a unique AGET chip for a slow control command.
 * @param asadIdx Index of the AsAd board to select (between 0 and 3).
 * @param agetIdx Index of the AGET chip on its AsAd board (between 0 and 3).
 */
void CoBoNode::selectAgetForSlowControl(size_t asadIdx, const size_t agetIdx)
{
	// Select AGET to write to
	const uint16_t mask = (1ul << (4*asadIdx+agetIdx));
	setAgetSlowControlMask(mask);
}
// --------------------------------------------------------------------------------------------------------------------
void CoBoNode::impl_prepare(CCfg::CConfig nodeCfg)
{
	try
	{
		// Remember which AsAd boards and AGET chips are active
		setAgetMask(nodeCfg);

		if (isCoBoActive())
		{
			NODE_INFO() << "Preparing: " << id().toString() << "'...";
			// Set CoBo's initDone bit to 0 so that CoBo board sends 'not configured' status to Mutant.
			setInitDone(false);
			// Check hardware version and firmware versions
			checkReleaseNumbers();
			// Initialize firmware pipeline
			pipeInit();
			// Initialize LEMO connectors
			setLemoModes(nodeCfg);
			// Initialize trigger mode
			triggerInit();
			// Initialize circular buffer
			circularBufferInit();
			// Initialize AsAd registers
			asadInit(nodeCfg, 10);
			// Initialize CoBo clock frequencies
			clockInit(nodeCfg);
			// Set manual SCR data delay
			setReadDataDelay(nodeCfg);
			// Set AsAd AGET input protections
			asadAgetInputManagerInit(nodeCfg);
			// Initialize AGET
			agetInit(nodeCfg);
			// Initialize ADC parameters
			asadAdcInit(DDR);
			// Initialize AsAd DAC registers
			asadDacInit();
			// Initialize AsAd auto-protection
			asadAlertsInit(nodeCfg);
			// Subscribe to alarm notifications
			subscribeToCoBoAlarms(nodeCfg);
			asadMonitoringRead();
			// Initialize value of AGET global threshold (why now?)
			setGlobalRegThresholds(UINT32_C(2));
			// Compute automatic SCR data delay
			calibrateAutoReadDataDelay(nodeCfg);
			// Select automatic or manual SCR data delay
			selectReadDataDelay(nodeCfg);
		}
	}
	catch (const CCfg::Exception & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Powers off AsAd boards, un-subscribes from alarms.
 */
void CoBoNode::impl_unprepare()
{
	try
	{
		// Check monitoring status of active AsAd boards (it will no longer be possible after they are powered off).
		uint32_t asadEnable = ctrl()->readRegister("asadEnable");
		for (size_t asadId = 0; asadId < 4u; asadId++)
		{
			if (not ::utl::BitFieldHelper<uint32_t>::getBit(asadEnable, asadId))
				continue;
			selectAsAdForRead(asadId);
			const uint16_t status = asad()->readRegister("interruptStatus1");
			if (0 != status)
			{
				NODE_WARN() << "AsAd board no. " << asadId << " is in auto-protection mode (status "
					<< std::hex << std::showbase << status << std::dec << " e.g. "
					<< get::asad::utl::stringifyMonitoringAlertStatus(status) << ')';
			}
		}

		// Power off AsAd boards
		ctrl()->writeRegister("asadEnable", 0);

		// Un-subscribe from alarms
		const mt::AlarmServicePrx& alarmSvce = ctrlPrx()->getAlarmService();
		alarmSvce->unsubscribe("ECC");
	}
	catch (const CmdException & e)
	{
		NODE_WARN() << "Error while unpreparing CoBo board (" << e.reason << ")";
	}
	catch (const Ice::LocalException & e)
	{
		NODE_WARN() << "Error while unpreparing CoBo board (" << e.ice_name() << ")";
	}
	catch (const Ice::UserException & e)
	{
		NODE_WARN() << "Error while unpreparing CoBo board (" << e.ice_name() << ")";
	}
}
// --------------------------------------------------------------------------------------------------------------------
void CoBoNode::impl_configure(CCfg::CConfig nodeCfg)
{
	try
	{
		// Remember which AsAd boards and AGET chips are active
		setAgetMask(nodeCfg);
		if (not isCoBoActive()) return;

		NODE_INFO() << "Configuring CoBo board '" << id().toString() << "'...";
		// Loop over each of the 4 AsAd boards
		for (size_t asadId=0; asadId < 4u; ++asadId)
		{
			// Skip inactive AsAd board
			if (isAsadActive(asadId))
			{
				// Loop over AGET chips
				switchToSlowControl();
				for (size_t agetId = 0; agetId < 4u; ++agetId)
				{
					if (not isAgetActive(asadId, agetId))
						continue;

					// Select AGET chip for slow control
					selectAgetForSlowControl(asadId, agetId);

					// Set AGET chip registers
					agetSetRegisters(nodeCfg, asadId, agetId);
					// Set CoBo registers for AGET chip 'readAlways' and 'readIfHit' masks
					agetSetReadIfHitMask(nodeCfg, asadId, agetId);
					agetSetReadAlwaysMask(nodeCfg, asadId, agetId);
					// Set CoBo thresholds for zero suppression mode
					setCoBoZeroSuppressionThresholds(nodeCfg, asadId, agetId);
				}
				switchToSlowControl();

				// Select AsAd board for slow control
				selectAsAdForSlowControl(asadId);

				// Set AsAd generator mode
				configureAsAdTcmDevice(nodeCfg, asadId);

				// Set AsAd inspection manager
				setAsAdInspectionManager(nodeCfg, asadId);
			}
		}

		// Reset AsAd and AGET masks
		setAsAdSlowControlMask(asadMask_);
		setAgetSlowControlMask(agetMask_);

		// Enable/disable readout of AGET FPN channels
		coboEnableFPNRead(nodeCfg);
		// Enable/disable hit register modification
		setWriteHitRegisterMask(nodeCfg);
		// Enable/disable CoBo test acquisition mode
		setTestMode(nodeCfg);
		// Set parameters of multiplicity trigger mode
		setMultiplicityMode(nodeCfg);
		setMultiplicityThreshold2p(nodeCfg);
		// Set multiplicity trigger moving average/deviation subtraction
		setCoBoMultSubtractMode(nodeCfg);
		// Select which chips to ignore in the multiplicity trigger
		setCoBoSuppressMultiplicity(nodeCfg);
		// Enable/disable 2p mode in CoBo
		enableMem2pMode(nodeCfg);
		// Set data source ID
		setDataSourceId(nodeCfg);
		// Set CoBo ID
		setCoBoId(nodeCfg);
		// Set trigger period, delays, ...
		setTriggerPeriod(nodeCfg);
		setTriggerDelay(nodeCfg);
		setTriggerDeadTime2p(nodeCfg);
		setTriggerDelay2p(nodeCfg);
		setTriggerTime2p(nodeCfg);
		setTriggerTimeOut2p(nodeCfg);
		// Set circular buffer parameters
		setCircularBuffers(nodeCfg);
		// Select trigger mode
		setTriggerModes(nodeCfg);
		setScwMultDelay(nodeCfg);
		setScwScrDelay(nodeCfg);
		// Set CoBo readout mode
		setReadoutMode(nodeCfg);
		// Set CoBo zero suppression mode
		setCoBoZeroSuppressionMode(nodeCfg);
		// Set number of AGET time buckets to read
		setReadoutDepth(nodeCfg);
		// Set readout offset in data frame header
		setReadoutOffset(nodeCfg);
		// Setup AsAd periodic pulser
		setAsAdPeriodicPulser(nodeCfg);

		// Initialize 2p mode
		aget2pInitNow();

		// Remember whether to reset timestamp and event counter during start
		try
		{
			useExternalTimeReset_ = nodeCfg("Module")("useExternalTimeReset").getBoolValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			useExternalTimeReset_ = false;
		}
	}
	catch (const CCfg::Exception & e)
	{
		throw CmdException(e.what());
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Read temperature of all active AsAd boards
 */
void CoBoNode::readAsAdTemperatures()
{

	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		if (not isAsadActive(asadIdx)) continue;
		const int16_t Ti = readAsAdIntTemp(asadIdx);
		const int16_t Te = readAsAdExtTemp(asadIdx);
		NODE_INFO() << "Read temperatures of AsAd board no. " << asadIdx << " of CoBo board '" << id().toString() << "': Tint = " << Ti << " deg C, Text = " << Te << " deg C";
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Performs pre-run actions.
 */
void CoBoNode::runInit()
{
	if (not isCoBoActive()) return;

	// Read temperature of all active AsAd boards
	readAsAdTemperatures();

	// Reset timestamp and event counters
	if (not useExternalTimeReset_)
	{
		NODE_DEBUG() << "Resetting timestamp and event counters";
		ctrl()->writeRegister("resetTime", 0);
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Resets each read pointer to the beginning of its circular buffer.
 */
void CoBoNode::resetCoBoReadPointers()
{
	try
	{
		for (size_t asadIdx=0; asadIdx <4u; ++asadIdx)
		{
			std::string prefix = boost::lexical_cast< std::string >(asadIdx);
			ctrl()->writeRegister(std::string("readPtr_asad") + prefix, ctrl()->readRegister(std::string("memStart_asad") + prefix));
		}
	}
	catch (const CmdException &)
	{
		// Case of single AsAd firmware
		ctrl()->writeRegister(std::string("readPtr"), ctrl()->readRegister(std::string("memStart")));
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Starts transfer of data present in the embedded circular buffer.
 */
void CoBoNode::impl_start()
{
	if (isCoBoActive())
	{
		hwNode().regProc("ctrl")
			.write("pipeCtrl", "firmwareEnable", 0) // Make sure firmware is not enabled
			.write("pipeCtrl", "intEnable", 0) // Make sure firmware does not generate interrupts
			// Resets firmware FIFO and state machine
			.write("pipeCtrl", "reset", 1) // Resets write pointers to memStart values
			.write("pipeCtrl", "backpressure", 1); // obsolete ?

		resetCoBoReadPointers(); // Resets read pointers to value of write pointers

		// Start data sender and interrupt monitor
		if (not hwNode().dataSenderCtrlPrx)
		{
			throw CmdException("Cannot start run without a valid data sender!");
		}
		hwNode().dataSenderCtrlPrx->daqStart();

		hwNode().regProc("ctrl")
			// Enable firmware for acquisition
			.write("pipeCtrl", "backpressure", 0) // obsolete ?
			.write("pipeCtrl", "reset", 0) // firmware is now in post-reset state
			// Firmware can now issue interrupts
			.write("pipeCtrl", "intEnable", 1)
			.write("pipeCtrl", "firmwareEnable", 1); // go!

		switchToFastControl();

		postStart();
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Actions to be run after the start of the run.
 */
void CoBoNode::postStart()
{
	try
	{
		if (configData_.asadPulserPeriodicModeEnabled)
		{
			NODE_INFO() << "Starting periodic AsAd pulses...";
			asadPulserMgr()->startPeriodicPulser();
		}
	}
	catch (mdaq::utl::CmdException & e)
	{
		NODE_ERROR() << e.reason;
	}
	catch (std::exception & e)
	{
		NODE_ERROR() << e.what();
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Actions to be run before stopping the run.
 */
void CoBoNode::preStop()
{
	try
	{
		// Stop periodic AsAd pulser
		asadPulserMgr()->stopPeriodicPulser();

		// Read temperature of all active AsAd boards
		readAsAdTemperatures();
	}
	catch (mdaq::utl::CmdException & e)
	{
		NODE_ERROR() << e.reason;
	}
	catch (std::exception & e)
	{
		NODE_ERROR() << e.what();
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Stops transfer of data present in the embedded circular buffer.
 */
void CoBoNode::impl_stop()
{
	if (not isCoBoActive()) return;

	try
	{
		preStop();

		// Disable firmware acquisition
		ctrl()->writeField("pipeCtrl", "firmwareEnable", 0);

		// Stop interrupts
		pipeInterruptEnable(false);

		// Stop listening to interrupts and flush data
		if (hwNode().dataSenderCtrlPrx)
			hwNode().dataSenderCtrlPrx->daqStop();
	}
	// Protect against call on rebooted CoBo board ('ctrl' device no longer exist).
	catch (const Ice::LocalException & e)
	{
		NODE_WARN() << "Error stopping CoBo board (" << e.ice_name() << ')';
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Connects to data routers.
 */
void CoBoNode::daqConnect()
{
	// Connect sender(s) to router(s)
	Node::daqConnect();

	// Select which circular buffers to read
	setCircularBuffersEnabled(asadMask_);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Reads the 'hardwareVersion' CoBo register.
 *
 * This 'hardwareVersion' register (offset 0x1D4) was introduced with firmware 20131107.
 * It has two valid values: 1= Reduced CoBo, >=2= CoBo v0.2.
 *
 * @return Value of hardware version register or 1 if not found or invalid.
 */
uint32_t CoBoNode::readHardwareVersion()
{
	Ice::Long versionId = 1;
	try
	{
		hwNode().regProc("ctrl").read("hardwareVersion", versionId);
		// Firmware version < 20131107 did not set 'hardwareVersion' register. Assume R-CoBo if register value is invalid.
		if (0xffffffffLL == versionId)
		{
			versionId = 1;
		}
	}
	catch (const CmdException &)
	{
		// Assume R-CoBo if register has not been not described
	}
	return versionId;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Checks firmware release number and release date.
 *
 * They are provided by CoBo registers 'firmwareDate' (offset 0x1D8) and 'firmwareRelease' (offset 0x1DC).
 * firmwareDate is a hexadecimal value displaying the build date of the firmware in the format 0xYYYYMMDD (ie. firmware from Nov 7th, 2013 will be value 0x20131107).
 * firmwareRelease stores the sequential release number of the firmware.  The firmware build on Nov 7th, 2013 will be version 1 since it is the first version with this register.
 *
 * @return Firmware release number, or 0 if register has not been found.
 */
uint32_t CoBoNode::checkFirmwareRelease()
{
	std::ostringstream oss;
	oss << "Checking firmware release... ";
	uint32_t firmwareRelease = 0;
	try
	{
		firmwareRelease = ctrl()->readRegister("firmwareRelease");
		const uint32_t firmwareDate = ctrl()->readRegister("firmwareDate");
		oss << firmwareRelease << " (" << std::hex << firmwareDate << ')';
		NODE_DEBUG() << oss.str();
	}
	catch (const CmdException &)
	{
		oss << " not found";
		NODE_WARN() << oss.str();
	}
	return firmwareRelease;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Checks hardware version and firmware release.
 */
void CoBoNode::checkReleaseNumbers()
{
	// Check firmware release
	checkFirmwareRelease();

	// Check hardware version
	std::ostringstream oss;
	oss << "Checking hardware version... ";
	const uint32_t hardwareVersion = readHardwareVersion();
	oss << hardwareVersion << " (" << ((2 <= hardwareVersion)?"CoBo >= v2.0":"R-CoBo") << ')';
	NODE_DEBUG() << oss.str();
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Initializes the 'pipeCtrl' register of the current CoBo node.
 * This registers controls the firmware pipeline.
 */
void CoBoNode::pipeInit()
{
	NODE_DEBUG() << "Initializing CoBo firmware pipeline";
	hwNode().regProc("ctrl")
		.write("pipeCtrl", "multiplicityReadOut", 0)   // Disable creation of frames with multiplicity data (depreciated, specific to 1 AGET).
		.write("pipeCtrl", "modeEnable_2p", 0)         // Disable 2p mode.
		.write("pipeCtrl", "bitSlipCalibrate", 0)      // Disable CoBo deserializer calibration sequence.
		.write("pipeCtrl", "partialReadOut", 0)        // Disable partial readout mode.
		.write("pipeCtrl", "testMode", 0)              // Disable test mode (where frames contain incremental counter instead of actual data).
		.write("pipeCtrl", "zeroSuppressionEnable", 0) // Disable zero suppression mode.
		.write("pipeCtrl", "firmwareEnable", 0)        // Disable the firmware state machine.
		.write("pipeCtrl", "hitRegisterOverride", 0)   // Depreciated.
		// Corresponds to script sequencer.ecc
		.write("pipeCtrl", "SCW", 0)
		.write("pipeCtrl", "SCR", 0)
		.write("pipeCtrl", "ScrScwOverride", 1)        // Allow software override of the SCW and SCR bit fields.
		.write("pipeCtrl", "reset", 1)                 // Reset all the CoBo logic.
		.write("pipeCtrl", "reset", 0)
		.write("pipeCtrl", "backpressure", 0);         // Depreciated (specific to 1 AGET).

	pipeInterruptEnable(false);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Allows or disables interrupts.
 */
void CoBoNode::pipeInterruptEnable(bool enable)
{
	hwNode().regProc("ctrl")
		.write("pipeCtrl", "intEnable", enable?1:0)
		.sleep_ms(500);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Initializes the trigger.
 */
void CoBoNode::triggerInit()
{
	NODE_DEBUG() << "Initializing CoBo trigger core";

	// Set manual delay between rising edge of SCR and first data on CoBo
	// 15 corresponds to a 1 meter cable
	// 16 corresponds to a 3 meter cable (not sure)
	for (size_t agetId=0; agetId < 4u; agetId++)
	{
		ctrl()->writeField("readDataDelay", "AGET" + boost::lexical_cast< std::string >(agetId), 14);
	}

	// Delay (units of 10 ns) between event causing the trigger and actual freeze of the AGET memory.
	ctrl()->writeRegister("triggerDelay", 0); // delay before SCW goes to 0, in multiples of 10ns

	// Delay (units of 10 ns) between event causing the second trigger and actual freeze of the AGET second memory.
	ctrl()->writeRegister("triggerDelay_2p", 0);

	// Delay (units of 10 ns) after the first trigger after which the second trigger can be issued.
	ctrl()->writeRegister("triggerWindowStart_2p", 0);

	// Maximum duration (units of 10 ns) to wait for the second trigger. After this delay, first trigger is canceled.
	ctrl()->writeRegister("triggerWindowEnd_2p", 100000); // 1ms

	// Multiplicity threshold for the second trigger.
	ctrl()->writeRegister("multThreshold_2p", 3000);

	// Period of the second trigger if periodic, in units of 10 ns.
	ctrl()->writeRegister("triggerTime_10ns_2p", 1000); // 10us

	// Multiplicity threshold (total of the 4 AsAd's).
	ctrl()->writeRegister("multThreshold", 3000);

	// Period of the trigger, if periodic, in units of 10 ns.
	ctrl()->writeRegister("triggerPeriod_10ns", 500000000); // 5s

	// Trigger mode
	ctrl()->writeRegister("triggerMask", 0); // no trigger mode selected

	// Number of samples to sum in the multiplicity sliding window for the multiplicity mode.
	ctrl()->writeRegister("multWindSize", 16);

	// Number of time buckets to readout.
	ctrl()->writeRegister("readDepth", 512); // is set by default to the reading of all time buckets

	// Channels to always read in partial readout mode
	ctrl()->writeRegister("readAlwaysMask1_4", 0);
	ctrl()->writeRegister("readAlwaysMask5_36", 0);
	ctrl()->writeRegister("readAlwaysMask37_68", 0);

	// Channels to read only if hit in partial readout mode.
	ctrl()->writeRegister("readIfHitMask1_4", 0xF);
	ctrl()->writeRegister("readIfHitMask5_36", 0xFFFFFFFF);
	ctrl()->writeRegister("readIfHitMask37_68", 0xFFFFFFFF);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Initializes the circular buffer.
 */
void CoBoNode::circularBufferInit()
{
	NODE_DEBUG() << "Initializing CoBo circular buffer(s)";

	// Memory addresses of the circular buffers of the 4 AsAd boards
	// Total space = 0x8000000 to 0x10000000 = 128 MiB
	uint64_t memStart_asad0 = UINT64_C(0x8000000);
	uint64_t bufferSize_B = UINT64_C(0x2000000); // 32 MiB
	try
	{
		for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
		{
			uint64_t memStart = memStart_asad0 + asadIdx*bufferSize_B;
			// Start and end
			ctrl()->writeRegister(std::string("memStart_asad") + boost::lexical_cast< std::string >(asadIdx), memStart);
			ctrl()->writeRegister(std::string("memEnd_asad") + boost::lexical_cast< std::string >(asadIdx), memStart + bufferSize_B);
			// Read pointer
			ctrl()->writeRegister(std::string("readPtr_asad") + boost::lexical_cast< std::string >(asadIdx), memStart);
		}
	}
	// Case of old R-CoBo hardware description
	catch (const CmdException &)
	{
		// Start and end
		ctrl()->writeRegister(std::string("memStart") , memStart_asad0);
		ctrl()->writeRegister(std::string("memEnd"), memStart_asad0 + 4*bufferSize_B);
		// Read pointer
		ctrl()->writeRegister(std::string("readPtr"), memStart_asad0);
	}

	// Number of events per CPU interrupt
	try
	{
		ctrl()->writeRegister("interruptRate", 0x1);
	}
	catch (const CmdException &) {}

	// Threshold (Bytes) indicating the amount of free space below which the firmware stops writing in the circular buffer.
	// Should be greater than the maximum frame size.
	ctrl()->writeRegister("aFullLevel", 0x01000000);

	// Level (Bytes) below which the firmware resumes writing data to the circular buffer (see 'aFullLevel').
	ctrl()->writeRegister("aEmptyLevel", 0x01000000);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets all CoBo LEDs off.
 */
void CoBoNode::clearAllLeds()
{
	try
	{
		ledManager()->setLEDs(false);
	}
	catch (const Ice::OperationNotExistException & e)
	{
		NODE_WARN() << e.what();
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Modifies the state of a given CoBo LED.
 */
void CoBoNode::modifyLed(const cobo::LedType & type, size_t asadId, const cobo::LedState & state)
{
	try
	{
		ledManager()->modifyLED(type, asadId, state);
	}
	catch (const Ice::OperationNotExistException & e)
	{
		NODE_WARN() << e.what();
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Powers ON AsAd board(s) and sets the configuration bits
 * of the AsAd slow controller and of AsAd to their initial values.
 */
void CoBoNode::asadInit(CCfg::CConfig nodeCfg, size_t delay_ms)
{
	NODE_DEBUG() << "Initializing AsAd board(s)";

	// Switch all LED's off
	clearAllLeds();

	// Resets the AsAd slow controller on CoBo
	ctrl()->writeField("asadController", "reset", 1); // If pulsed, resets error flag, busy flag and checksum error
	ctrl()->writeField("asadController", "manualCal", 0); // Value of the manual delay (units of 10 ns)
	ctrl()->writeField("asadController", "calibration", 0); // If set, calibration is automatic (1 bit per AsAd)
	ctrl()->writeField("asadController", "reset", 0);

	// Powers ON AsAd boards
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		// Check AsAd monitoring device alert status register
		if (isAsAdConnected(asadIdx))
		{
			const uint16_t mask = (1ul << asadIdx);
			setAsAdSlowControlMask(mask);
			selectAsAdForRead(asadIdx);
			const uint16_t asadAlertStatus = asad()->readRegister("interruptStatus1");
			if (0 != asadAlertStatus)
			{
				NODE_WARN() << "AsAd board no. " << asadIdx << " of CoBo board '" << id().toString() << "' was in auto-protection mode: status="
						<< std::hex << std::showbase << asadAlertStatus << std::dec
						<< " (" << asad::utl::stringifyMonitoringAlertStatus(asadAlertStatus) << ')';
			}
		}

		std::string fieldName = std::string("powerON") + boost::lexical_cast< std::string >(asadIdx);
		// Power OFF board
		ctrl()->writeField("asadEnable", fieldName, 0);
		sleep_ms(delay_ms);
		if (isAsadActive(asadIdx))
		{
			// Check board is connected
			if (not isAsAdConnected(asadIdx))
			{
				modifyLed(cobo::LedB, asadIdx, cobo::LedSlowPulse);
				throw mdaq::utl::CmdException(std::string("AsAd board no. ") + boost::lexical_cast< std::string >(asadIdx) + std::string(" is not connected to CoBo!"));
			}
			// Check whether to keep going in case of problem
			bool abortIfNoPowerSupply = true;
			try
			{
				abortIfNoPowerSupply = asadConfig(nodeCfg, asadIdx)("Control")("checkPowerSupply").getBoolValue();
			}
			catch (const CCfg::CConfig::NotFound & e)
			{
				NODE_WARN() << e.what();
			}
			// Check board is power supplied
			if (not isAsAdPowerSupplied(asadIdx))
			{
				std::ostringstream errorMsg;
				errorMsg << "AsAd board no. " << asadIdx << " of CoBo board '" << id().toString() << "' seems disconnected from its power supply!";

				if (abortIfNoPowerSupply)
				{
					modifyLed(cobo::LedB, asadIdx, cobo::LedSlowPulse);
					throw mdaq::utl::CmdException(errorMsg.str());
				}
				else
				{
					NODE_ERROR() << errorMsg.str();
				}
			}
			// Power ON board
			NODE_DEBUG() << "Powering on AsAd board no. " << asadIdx;
			ctrl()->writeField("asadEnable", fieldName, 1);
			sleep_ms(delay_ms);
			// Check board status
			bool plg, alarm, pwu;
			asadPowerStatus(asadIdx, plg, alarm, pwu);
			if (alarm)
			{
				std::ostringstream msg;
				msg << "AsAd board no. " << asadIdx << " of CoBo board '" << id().toString() << "' was not successfully powered on"
						" (PLG=" << (not plg) << ", AL=" << ", PWU=" << pwu << "!";
				modifyLed(cobo::LedB, asadIdx, cobo::LedSlowPulse);
				throw mdaq::utl::CmdException(msg.str());
			}

			// Switch 'P' LED on to indicate AsAd board is successfully connected and powered.
			modifyLed(cobo::LedP, asadIdx, cobo::LedOn);
		}
	}

	// Set AsAd slow control mask
	setAsAdSlowControlMask(asadMask_);

	// Calibrates the AsAd slow control
	asadCal();

	// Configure and start AsAd monitoring (without auto-protection)
	asadMonitoringInit(nodeCfg);

	// Check AsAd VDD voltages
	asadMonitoringCheckVDD(nodeCfg);

	// Set AGET input protections
	asadSetAgetInputProtections();
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Initializes AsAd DAC device registers.
 * It seems the value of the DAC1 ("DATA") register is not always its official power-on value of 0.
 * @see GET-AS-002-0010.
 */
void CoBoNode::asadDacInit()
{
	asad()->writeRegister("DAC0", 0x81);
	asad()->writeRegister("DAC1", 0);
	// Register DAC2 (VERSION) is read only.
	asad()->writeRegister("DAC3", 0);
	asad()->writeRegister("DAC4", 0);
	asad()->writeRegister("DAC5", 0);
	asad()->writeRegister("DAC6", 0); // 0x38 according to GET-AS-002-0010
	asad()->writeRegister("DAC7", 0);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Calibrates the AsAd slow control.
 */
void CoBoNode::asadCal()
{
	NODE_DEBUG() << "Calibrating AsAd slow control ...";

	setAsAdSlowControlMask(asadMask_);

	// Synchronizes the slow control lines
	hwNode().regProc("ctrl")
		.write("asadController", "reset", 1) 			// Reset error flag, busy flag and checksum error.
		.write("asadController", "reset", 0)
		.write("asadController", "manualCal", 0x0)		// Value of the manual delay (units of 10 ns). Specific to Reduced CoBo.
		.write("asadController", "calibration", asadMask_)	// If set, calibration is automatic.
		.write("asad_sc_start", 0x1)				// Starts write/read operation on AsAd register.
		.write("asadController", "calibration", 0);

	// Check calibration of AsAd slow control was successful:
	// Manufacturer ID of monitoring device must be 0x41
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		// Skip unused board
		if (not isAsadActive(asadIdx)) continue;

		// Select AsAd board
		selectAsAdForSlowControl(asadIdx);
		// Read register
		uint64_t manufacturerID = asad()->readRegister("monitorID");

		// Check manufacturer ID value is correct
		if (manufacturerID != UINT64_C(0x41))
		{
			std::ostringstream msg;
			msg << "Calibration of slow control of AsAd no. " << asadIdx << " failed: expected 0x41, read "
					<< std::hex << std::showbase << manufacturerID << '!';
			throw mdaq::utl::CmdException(msg.str());
		}
		else
		{
			// Set user chosen bits of board identifier
			const uint8_t expectedUserID = 1 + asadIdx;
			asadSetBoardUserID(expectedUserID);

			// Read board identifier
			const uint32_t fullID = asadGetIdentifier();
			uint16_t actualUserID = utl::BitFieldHelper< uint32_t >::getField(fullID, 24u, 8u);
			uint32_t boardID = fullID & 0xFFFFFF;
			NODE_DEBUG() << "AsAd board no. " << asadIdx << " has identifier " << std::hex << boardID << std::dec
					<< " and user label #" << actualUserID;

			// Check user chosen bits
			if (actualUserID != expectedUserID)
			{
				std::ostringstream msg;
				msg << "Calibration of slow control of AsAd no. " << asadIdx << " failed: expected and actual user labels of AsAd no. " << asadIdx << " do not match!";
				throw mdaq::utl::CmdException(msg.str());
			}
		}
	}
	setAsAdSlowControlMask(asadMask_);
	NODE_DEBUG() << " ... AsAd slow control calibration done";
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Initializes writing CKW and reading CKR clock frequencies.
 */
void CoBoNode::clockInit(CCfg::CConfig nodeCfg)
{
	NODE_DEBUG() << "Initializing clock frequencies";
	setAsAdSlowControlMask(asadMask_);

	// CoBo CKW
	const float freqCKW_MHz = nodeCfg("Module")("writingClockFrequency").getRealValue();
	bool pllConfigEnabled = true;
	get::cobo::PllRegisterMap regs;
	try
	{
		pllConfigEnabled = not nodeCfg("Module")("skipPLLConfig").getBoolValue();
		// Check advanced PLL parameters
		for (size_t regIndex=10; regIndex <= 13; ++regIndex)
		{
			regs[regIndex] = rc::CoBoConfigParser::getPllRegValue(nodeCfg, regIndex);
		}
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		if (readHardwareVersion() >= 2)
		{
			NODE_WARN() << e.what(); // Debugging parameter added on April 22th, 2014
		}
	}
	setCoBoWritingClockFrequency(freqCKW_MHz, pllConfigEnabled, regs);

	// AsAd LCM1
	asadSamplingClockInit(nodeCfg, freqCKW_MHz);

	// CoBo CKR
	const float freqCKR_Mhz = nodeCfg("Module")("readingClockFrequency").getRealValue();
	setCoBoReadingClockFrequency(freqCKR_Mhz);

	// AsAd LCM2
	asadReadingClockInit(nodeCfg, freqCKR_Mhz);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets function of the 4 LEMO connectors on the CoBo front panel.
 */
void CoBoNode::setLemoModes(CCfg::CConfig nodeCfg)
{
	if (readHardwareVersion() < 2u) return;

	NODE_DEBUG() << "Configuring LEMO connectors";
	try
	{
		for (size_t lemoIndex=0; lemoIndex<4u; lemoIndex++)
		{
			// Check desired function
			uint8_t lemoMode = getLemoConnectorMode(nodeCfg, lemoIndex);
			// Set register
			const std::string fieldName = std::string("L") + boost::lexical_cast< std::string >(lemoIndex);
			ctrl()->writeField("lemoMux", fieldName, lemoMode);
		}
	}
	catch (const mdaq::utl::CmdException & e)
	{
		NODE_WARN() << e.what();
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the @em initDone bit of CoBo's @em mutantConfig Mutant configuration register to 0 or 1.
 *
 * This bit must be set to 0 while slow control configures CoBo and to 1 afterward.
 *
 * Email from N. U., September 23rd, 2014:
 * The initDone bit only affects the status the CoBo sends to Mutant.
 * It determines whether the CoBo sends a @em ready status or a <em>not configured</em> status.
 * It should always be set to 0 as the first step of configuring the CoBo and set to 1 as the last step of configuring the CoBo.  The bit has no effect on the function of the slow control; you can still use eccClient to change register values even during readout.  The initDone bit has no effect when there is no Mutant in the system.
 *
 * @param done Whether to set the 'initDone' bit to 0 (false) or 1 (true).
 */
void CoBoNode::setInitDone(bool done)
{
	if (not isCoBoActive()) done = false;
	try
	{
		ctrl()->writeField("mutantConfig", "initDone", done);
	}
	catch (const Ice::LocalException & e)
	{
		NODE_WARN() << "Could not set CoBo 'initDone' bit to " << done << " (" << e.ice_name() << ')';
	}
	catch (const mdaq::utl::CmdException & e)
	{
		NODE_WARN() << "Could not set CoBo 'initDone' bit to " << done << ". " << e.reason << ')';
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the @em mutantMode bit of CoBo's @em mutantConfig Mutant configuration register to 0 or 1.
 *
 * This bit must be set by slow control when the CoBo is connected to a Mutant.
 *
 * Email from N. U., September 23rd, 2014:
 * The mutantMode bit should always be set to 1 when using a Mutant module
 *  and it should always be 0 when there is no Mutant in the system.
 * The timing of setting the bit should not matter,
 *  but I recommend setting it to 1 at the beginning of configuration,
 *  just after you set the initDone bit to 0.
 */
void CoBoNode::setMutantMode(bool enabled)
{
	if (not isCoBoActive()) enabled = false;
	try
	{
		ctrl()->writeField("mutantConfig", "mutantMode", enabled);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		NODE_WARN() << e.what();
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the @em l2Trigger bit of CoBo's @em mutantConfig Mutant configuration register to 0 or 1.
 */
void CoBoNode::setLevel2Trigger(bool enabled)
{
	if (not isCoBoActive()) enabled = false;
	try
	{
		ctrl()->writeField("mutantConfig", "l2Trigger", enabled);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		NODE_WARN() << e.what();
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configures register 'ckwDivider'of CoBo 'ctrl' device.
 * The registers 'ckwDivider' and 'ckrDivider' are not used in the µTCA CoBo (see PLL device instead).
 * @param frequency_MHz Sampling clock frequency.
 */
void CoBoNode::setCoBoCkwDivider(const float& frequency_MHz)
{
	// Check value is a valid CKW frequency , compute CoBo frequency divider
	uint8_t ckwDivider = 0;
	if (frequency_MHz == 100.0)
	{
		ckwDivider = UINT8_C(0);
	}
	else if (frequency_MHz == 50.0)
	{
		ckwDivider = UINT8_C(1);
	}
	else if (frequency_MHz == 25.0)
	{
		ckwDivider = UINT8_C(2);
	}
	else if (frequency_MHz == 12.5)
	{
		ckwDivider = UINT8_C(3);
	}
	else if (frequency_MHz == 6.25)
	{
		ckwDivider = UINT8_C(4);
	}
	else if (frequency_MHz == 3.125)
	{
		ckwDivider = UINT8_C(5);
	}
	else
	{
		throw CmdException(boost::lexical_cast<std::string>(frequency_MHz) + " MHz is not a valid CKW frequency for R-CoBo! CKW must be of the form 100 MHz/2**N.");
	}

	// Configure register 'ckwDivider'of CoBo 'ctrl' device (R-CoBo only)
	// CKW = 100 MHz / 2**ckdDivider

	hwNode().regProc("ctrl").write("ckwDivider", "d", ckwDivider);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets CoBo writing (sampling) clock CKW frequency.
 *   - sets the AsAd 'LCM1' register (WCKn sampling clock frequency and phase adjustment) of currently selected AsAd boards for given frequency
 *   - sets the CoBo 'ckwDivider' register for the CKW frequency divider: CKW = 100 MHz / 2**ckwDivider (R-CoBo only).
 *   - configures M-CoBo PLL device for given frequency.
 *
 * @see GET-AS-002-0006
 *
 * Email from N. Usher (September 2014):
 * Using the AsAd to divide the clock might cause some problems.
 * The CoBo synchronizes SCW with CKW, and must pulse SCW for a specific number of clock cycles during initialization.
 * If you divide CKW in AsAd, the pulses will not be the correct length and might not be synchronized properly when they reach AGET.
 *
 * @param node CoBo node.
 * @param freqCKW_MHz CoBo CKW frequency [MHz].
 * @param pllConfigEnabled Debugging flag to disable configuration of the CoBo PLL device.
 * @param regs Map of non-default PLL register values.
 */
void CoBoNode::setCoBoWritingClockFrequency(const float & freqCKW_MHz, bool pllConfigEnabled, const cobo::PllRegisterMap & regs)
{
	NODE_DEBUG() << "Setting CKW clock frequency to " << freqCKW_MHz << " MHz";

	// Check whether node is an R-CoBo or an M-CoBo node: 1 = R-CoBo, 2 = CoBo v2.0
	const uint32_t hardwareVersion = readHardwareVersion();

	// For M-CoBo, configure the PLL device
	if (pllConfigEnabled and hardwareVersion >= 2)
	{
		CoBoPLLManager(*this).loadPLLConfig(freqCKW_MHz, regs);
	}
	// For R-CoBo, configure register 'ckwDivider'of CoBo 'ctrl' device
	else
	{
		setCoBoCkwDivider(freqCKW_MHz);
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets CoBo reading clock CKR frequency.
 *
 * Sets the CoBo 'ckrDivider' register for the CKR frequency divider: CKR = 100 MHz / 2**ckrDivider
 * The registers 'ckwDivider' and 'ckrDivider' are not used in the µTCA CoBo (see PLL device instead).
 *
 * @param frequency_MHz Frequency [MHz].
 */
void CoBoNode::setCoBoReadingClockFrequency(const float & frequency_MHz)
{
	// The register 'ckrDivider' is not used in the µTCA CoBo (see PLL device instead).
	if (readHardwareVersion() >= 2) return;

	NODE_DEBUG() << "Setting CKR clock frequency to " << frequency_MHz << " MHz";
	if (frequency_MHz == 25.0 or frequency_MHz == 12.5)
	{
		hwNode().regProc("ctrl").write("ckrDivider", "d", 2); //  25 MHz
	}
	else
	{
		std::ostringstream stream;
		stream << frequency_MHz << " MHz is not a valid reading clock (CKR) frequency. The CKR frequency must always be set to 25 MHz.";
		throw CmdException(stream.str());
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Adjusts the value of the LCM1 register of all active AsAd boards using the configuration file.
 * @param nodeCfg CoBo board configuration.
 * @param freqCKW_MHz Frequency of the CKW sampling clock provided by CoBo.
 */
void CoBoNode::asadSamplingClockInit(CCfg::CConfig nodeCfg, const float & freqCKW_MHz)
{
	for (size_t asadId = 0; asadId < 4; ++asadId)
	{
		if (isAsadActive(asadId))
		{
			// Compute LCM1 value
			float samplingClockDivider = 1; // WCKn = CKW
			try
			{
				samplingClockDivider = nodeCfg("AsAd", asadId)("Clocking")("samplingClockDivider").getRealValue();
			}
			catch (const CCfg::CConfig::NotFound & e)
			{
				NODE_WARN() << e.what();
			}
			uint64_t lcm1Value = asad::utl::buildLcm1RegisterValue(freqCKW_MHz, samplingClockDivider);
			// Patch LCM1 value
			const bool lcm1DebugMode = nodeCfg("AsAd", asadId)("Clocking")("Debug")("LCM1")("debugModeEnabled").getBoolValue(false);
			if (lcm1DebugMode)
			{
				::get::asad::utl::Lcm1Config lcm1Config;
				asad::utl::parseAsAdLcm1DebugConfig(nodeCfg("AsAd", asadId), lcm1Config);
				asad::utl::patchLcm1RegisterValue(lcm1Value, lcm1Config);
			}
			// Select AsAd
			selectAsAdForSlowControl(asadId);
			// Configure register LCM1 of AsAd Local Clock Manager device
			NODE_DEBUG() << "Setting AsAd no. " << asadId << " LCM1 register to "  << std::hex << std::showbase << lcm1Value << std::dec << std::noshowbase;
			asad()->writeRegister("LCM1", lcm1Value);
		}
	}
	setAsAdSlowControlMask(asadMask_);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Adjusts the value of the LCM2 register of all active AsAd boards using the configuration file.
 */
void CoBoNode::asadReadingClockInit(CCfg::CConfig nodeCfg, float freqCKR_Mhz)
{
	//NODE_DEBUG() << "Adjusting AsAd LCM2 configuration";
	for (size_t asadId = 0; asadId < 4; ++asadId)
	{
		if (isAsadActive(asadId))
		{
			// Compute LCM2 value
			int phaseShiftAdcClock_deg = nodeCfg("AsAd", asadId)("Clocking")("phaseShiftAdcClock").getIntValue();
			int delayAdcClock_ps = nodeCfg("AsAd", asadId)("Clocking")("delayAdcClock").getIntValue();
			uint64_t lcm2Value = asad::utl::buildLcm2RegisterValue(freqCKR_Mhz, phaseShiftAdcClock_deg, delayAdcClock_ps);
			// Patch LCM2 value
			const bool lcm2DebugMode = nodeCfg("AsAd", asadId)("Clocking")("Debug")("LCM2")("debugModeEnabled").getBoolValue(false);
			if (lcm2DebugMode)
			{
				::get::asad::utl::Lcm2Config lcm2Config;
				asad::utl::parseAsAdLcm2DebugConfig(nodeCfg("AsAd", asadId), lcm2Config);
				asad::utl::patchLcm2RegisterValue(lcm2Value, lcm2Config);
			}
			// Select AsAd
			selectAsAdForSlowControl(asadId);
			// Set LCM2 register
			NODE_DEBUG() << "Setting AsAd no. " << asadId << " LCM2 register to " << std::hex << std::showbase << lcm2Value << std::dec << std::noshowbase;
			asad()->writeRegister("LCM2", lcm2Value);
		}
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the manual read data delay.
 * Only has an effect when using a Reduced-CoBo board.
 * Manual delay registers have been disconnected on the micro-TCA CoBo.
 */
void CoBoNode::setReadDataDelay(CCfg::CConfig nodeCfg)
{
	int delay = nodeCfg("Module")("readDataDelay").getIntValue();
	NODE_DEBUG() << "Setting manual read data delay to " << delay;
	ctrl()->writeField("readDataDelay", "AGET0", delay);
	ctrl()->writeField("readDataDelay", "AGET1", delay);
	ctrl()->writeField("readDataDelay", "AGET2", delay);
	ctrl()->writeField("readDataDelay", "AGET3", delay);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Initializes AGET chip and its slow control.
 */
void CoBoNode::agetInit(CCfg::CConfig nodeCfg)
{
	NODE_DEBUG() << "Initializing AGET registers";

	// Reset AGet controller error and busy flags
	ctrl()->writeField("agetController", "reset", 0);
	ctrl()->writeField("agetController", "reset", 1);
	ctrl()->writeField("agetController", "reset", 0);

	// Set mask for broadcast of AGET slow control commands
	ctrl()->writeField("aget_sc_configuration", "SC_En", agetMask_);
	ctrl()->readField("aget_sc_configuration", "SC_En");

	// Set mask for broadcast of AsAd slow control commands
	setAsAdSlowControlMask(asadMask_);

	// Reset the AsAd state machine for the AGET slow control
	asad()->writeRegister("resetFast", 0); // Reset internal state machine of Fast Control process

	agetCal(nodeCfg);

	aget()->writeField("reg1", "DAC", 1);
	ctrl()->writeRegister("initSCAReadOutAfterPowerOnNow", 1);

	//DEBUG_END() << "agetInit()";
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets which AGET should be used for cable delay calibrations.
 *
 * The 'calibrationAget' CoBo register (offset 0x1E0) sets which AGET should be used for cable delay calibrations.
 * It has 4 fields, one for each AsAd.  Each field is two bits.
 * This register must be set before the calibration sequence is started.
 * This feature was added on Nov 7, 2013 (at the same time as the 'firmwareDate' register).
 * This setting should only be modified on AsAds that do not have an AGET installed in position 0.
 *
 * @param agetMask Mask of AGET chips that are available for calibration.
 * @param nodeCfg configuration to lookup which chips the user wants to use for calibration.
 */
void CoBoNode::coboSetCalibrationAgetChips(uint16_t agetMask, CCfg::CConfig nodeCfg)
{
	NODE_DEBUG() << "Selecting calibration AGET chips (mask=" << std::hex << std::showbase << agetMask << std::dec << ')';
	try
	{
		// Check if users wants to use specific chips for calibration
		const uint8_t userCalibMask = getCalibrationChips(nodeCfg);
		// Loop over AsAd boards
		for (size_t asadIdx=0; asadIdx < 4; ++asadIdx)
		{
			size_t calibrationChipId = 0;
			if (isAsadActive(asadIdx))
			{
				// Find first active chip of the 4 chips on this AsAd board
				for (calibrationChipId=0; calibrationChipId < 4; ++calibrationChipId)
				{
					if (utl::BitFieldHelper< uint16_t >::getBit(agetMask, asadIdx*4+calibrationChipId))
						break;
				}
				// Check which chip the user wants to use
				size_t userCalibrationChipId = utl::BitFieldHelper< uint8_t >::getField(userCalibMask, asadIdx*2, 2);
				if (userCalibrationChipId > 0)
				{
					calibrationChipId = userCalibrationChipId;
				}
				// Handle illegal case (no active chip on the board or erroneous configuration file)
				if (calibrationChipId > 3)
				{
					calibrationChipId = 0;
				}
			}
			if (calibrationChipId > 0)
			{
				NODE_INFO() << "Selecting AGET chip no. " << asadIdx*4+calibrationChipId << " for AsAd board no. " << asadIdx << " serial control delay calibration";
			}
			// Set AGET chip for serial control delay calibration
			std::string fieldName = std::string("AsAd") + boost::lexical_cast< std::string >(asadIdx);
			//NODE_DEBUG() << "Setting calibrationAget." << fieldName << " to " << calibrationChipId;
			ctrl()->writeField("calibrationAget", fieldName, calibrationChipId);
		}
	}
	catch (const CmdException & e)
	{
		// Case of old firmware/hardware description where the 'calibrationAget' register does not exist.
		NODE_WARN() << e.reason;
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Subscribes ECC to CoBo alarm notifications
 *  and enables CoBo monitoring of AsAd alarms.
 */
void CoBoNode::subscribeToCoBoAlarms(CCfg::CConfig & nodeCfg)
{
	try
		{
			// Check whether to monitor AsAd alarms
			bool checkAsAd = true;
			try
			{
				checkAsAd = nodeCfg("Alarms")("checkAsAd").getBoolValue();
			}
			catch (const CCfg::CConfig::NotFound & e)
			{
				NODE_WARN() << e.what();
			}
			// Enable monitoring of AsAd alarms
			NODE_DEBUG() << "Enabling periodic check of AsAd alarm bits by CoBo";
			asadAlarmMonitor()->setAsAdAlarmMonitoringEnabled(checkAsAd);

			// Remove all subscribers
			const mt::AlarmServicePrx& alarmSvce = ctrlPrx()->getAlarmService();
			alarmSvce->reset();
			// Check whether to subscribe
			bool subscribe = true;
			try
			{
				subscribe = nodeCfg("Alarms")("subscribe").getBoolValue();
			}
			catch (const CCfg::CConfig::NotFound & e)
			{
				NODE_WARN() << e.what();
			}
			// Subscribe
			if (subscribe)
			{
				NODE_DEBUG() << "Subscribing to CoBo alarm notifications";
				if (backend)
				{
					backend->selectHwNodeByEndpoint(hwNode().endpoint);
					alarmSvce->subscribe("ECC", backend->alarmLogger());
				}
			}
		}
		catch (CmdException & e)
		{
			NODE_WARN() << e.reason;
		}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Calibrates the AGET slow control
 *
 * By default, during the calibration sequence, the firmware measures the time delay for data received from AGET 0.
 */
void CoBoNode::agetCal(CCfg::CConfig nodeCfg)
{
	NODE_DEBUG() << "Calibrating AGET slow control";

	// Create local mask for easier debugging
	const uint16_t mask = agetMask_;

	// Select which AGET to use for serial control delay calibration
	coboSetCalibrationAgetChips(agetMask_, nodeCfg);

	adcReset();
	adcTest(DDR);

	ctrl()->writeField("ckpFastDivider", "d", 0);
	ctrl()->writeField("agetController", "reset", 1);
	ctrl()->writeField("agetController", "reset", 0);

	ctrl()->writeField("aget_sc_configuration", "SC_En", mask);
	ctrl()->writeField("aget_sc_configuration", "pulseEnable", mask);

	ctrl()->writeField("agetController", "hitRegister", 0);
	ctrl()->writeField("agetController", "calibration", 0);
	//ctrl()->writeRegister("aget_sc_rw", 1); // < from agetCal.ecc
	sleep_ms(200);

	aget()->readRegister("reg5");
	ctrl()->writeField("agetController", "reset", 1);
	ctrl()->writeField("agetController", "reset", 0);

	ctrl()->writeField("aget_sc_configuration", "pulseEnable", mask);

	ctrl()->writeField("agetController", "hitRegister", 0);
	ctrl()->writeField("agetController", "calibration",1);

	aget()->readRegister("reg5");
	ctrl()->writeField("agetController", "reset", 1);
	ctrl()->writeField("agetController", "reset", 0);
	ctrl()->writeField("agetController", "calibration", 0);

	ctrl()->writeField("aget_sc_configuration", "pulseEnable", mask);

	ctrl()->writeField("agetController", "hitRegister", 1);
	aget()->writeRegister("hit", UINT64_C(0xA000000000000005));
	ctrl()->writeField("agetController", "calibration", 1);

	ctrl()->writeField("aget_sc_configuration", "pulseEnable", 0);

	ctrl()->writeField("agetController", "hitRegister",1);
	aget()->writeRegister("hit", UINT64_C(0xA000000000000005));
	ctrl()->writeField("agetController", "calibration", 0);

	ctrl()->writeField("aget_sc_configuration", "pulseEnable", mask);
	ctrl()->writeField("agetController", "hitRegister", 0);

	aget()->readRegister("reg5");

	ctrl()->writeField("agetController", "reset", 1);
	ctrl()->writeField("agetController", "reset", 0);

	// Print AGET chip version numbers so that user can see if calibration of AGET slow control was successful
	printAgetChipVersionNumbers();
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Print version numbers of active AGET chips
 */
void CoBoNode::printAgetChipVersionNumbers(bool check)
{
	// Switch to slow control mode to be able to read AGET register 5
	switchToSlowControl();

	if (check)
	{
		NODE_DEBUG() << "Checking AGET chip version numbers are valid...";
	}

	std::ostringstream msg;
	msg << "AGET chip version numbers:";
	for (size_t asadId = 0; asadId < 4u; ++asadId)
	{
		for (size_t agetId=0; agetId < 4u; ++agetId)
		{
			if (utl::BitFieldHelper< uint16_t >::getBit(agetMask_, 4*asadId+agetId))
			{
				// Select AGET to read from
				ctrl()->writeField("agetChipSelect", "d", 4 * asadId + agetId);
				// Read register 5
				const uint16_t agetVersion = (uint16_t) aget()->readRegister("reg5");
				msg << ' ' << std::showbase << std::hex << agetVersion << std::dec;
				// Check version is known
				if (check and std::find(KNOWN_AGET_VERSIONS, KNOWN_AGET_VERSIONS+NUM_KNOWN_AGET_VERSIONS, agetVersion) == KNOWN_AGET_VERSIONS+NUM_KNOWN_AGET_VERSIONS)
				{
					msg.str("");
					msg << "Read unknown version '" << std::showbase << std::hex << agetVersion << std::dec << " for AGET chip no. " << agetId << " of AsAd board no. " << asadId;
					msg << " It means calibration of AGET slow control failed!";
					throw CmdException(msg.str());
				}
			}
			else
			{
				msg << " N/A";
			}
		}
	}
	NODE_INFO() << msg.str();
}
// --------------------------------------------------------------------------------------------------------------------
/**
 *  Triggers AsAd ADC software reset - resets all internal registers.
 */
void CoBoNode::adcReset()
{
	asad()->writeRegister("ADC0", 0x400);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the ADC to send a constant vector.
 */
void CoBoNode::adcTest(AdcMode adcMode)
{
	asad()->writeRegister("ADC2", 0x2A0);
	asad()->writeRegister("ADC3", 0x001);
	asad()->writeRegister("ADC4", 0x000);
	asad()->writeRegister("ADC5", (adcMode == SDR)?0x403:0x401);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Initializes AsAd ADC parameters.
 * See email from N. Usher on July 13th, 2013.
 * See Texas Instruments ADS6422 Data Sheet.
 * See email from J. Pibernat on October 16th, 2013.
 *
 * The ADC is configured through 8 write-only 11-bit registers (A to H, or ADC0 to ADC7).
 * Registers B, G and H (ADC1, ADC6 and ADC7) should remain at their initial value of 0x000.
 *
 *   -# Writing 0x400 to register A (ADC0) triggers a reset of all registers to the value 0x000 (including register A).
 *      Register A is also used to select which of the 4 ADC in the chip to switch off.
 *
 *   -# The purpose of register B (ADC1) is to reshape the ADC coding clock.
 *
 *   -# Le registre C (ADC2) va permettre de calibrer le processus de déserialisation.
 *      Comment cela ce passe? En fait, on ne met rien à l'entrée de l'ADC mais on force ses sorties à prendre certains états.
 *      Par exemple, en écrivant 0x0C3 on va sortir sur chacune des lignes qui sortent 6 bits une suite de 1 et de 0. Nathan va se servir de cela pour ajuster sa synchro.
 *      Si on écrit 0x2A0, alors on peut sortir le pattern qu'on veut sur les sorties.
 *
 *   -# The 12 bit custom pattern is defined through registers D (ADC3, lower 11 bits) and E (ADC4, most significant bit).
 *
 *   -# Register F (ADC5) specifies the output format of the ADC data.
 *      En écrivant 0x401 dans ce registre, on dit qu'une donnée codée sur 12 bits doit:
 *      - sortir sur deux voies qui distribuent chacune 6 bits en série
 *      - que sur une des deux voies on sort le MSB de la donnée (6 bits) et sur l'autre voie on sort le LSB (6 bits) avec pour chacun le bit de poids fort en premier.
 *
 *   -# Les registres G et H (ADC6 and ADC7) servent essentiellement à configurer les sorties d'ADC d'un point de vue "hardware"
 *      (pour bien comprendre ce que sont les "terminaisons" dans la datasheet, il faut savoir qu'en LVDS on ne travaille jamais en "tension").
 *       Une sortie LVDS différentielle est constituee de 2 lignes.
 *       Sur la première on pousse un courant i, sur la deuxième on tire le même courant.
 *       La conversion courant en tension se fait dans la charge qu'on appelle aussi la terminaison en appliquant la loi d'ohm tout bêtement U=RI.
 *       A l'entrée de CoBo il y a des R de 100 ohms (c'est le standard LVDS) et en sortie d'AsAd il y a les courants +i et -i qu'on envoie vers les terminaisons.
 *       L'ADC permet de "simuler" ces charges au cas où le récepteur  à distance (pour nous CoBo) ne serait pas connecté.
 *       Le réglage de ces charges se fait pas les 2 registres précédemment cités.
 *       Admettons qu'on ait fait ce réglage, si maintenant on branche CoBo on va se retrouver avec la charge "simulée" à la sortie de l'ADC en parallèle avec la 100 ohm de CoBo.
 *       Du coup, la résistance équivalente diminue.
 *       Si le courant reste constant ça signifie que la tension à l'entrée de CoBo sera plus faible. Pour maintenir cette tension constante sur une charge plus faible il faudrait augmenter i .
 *       C'est aussi ce que permettent de faire les registres G & H. Bref, on ne se sert pas de ces fonctions.
 */
void CoBoNode::asadAdcInit(AdcMode adcMode)
{
	NODE_DEBUG() << "Running AsAd ADC calibration sequence (" << (adcMode == DDR?"DDR":"SDR") << ')';

	// Set mask for broadcast of AsAd slow control commands
	setAsAdSlowControlMask(asadMask_);

	/* Reset ADC */
	adcReset();

	/* Calibrate of the ADC lines */
	adcCal(adcMode);

	/* Put the ADC in its default configuration */
	// Powers down ADC channels B, C, D
	//asad()->writeRegister("ADC0", 0x1C);
	// Straight binary data format, normal ADC operation
	asad()->writeRegister("ADC2", 0x200);
	// Output format of the ADC data: 2 Wire interface, bit clock selection
	asad()->writeRegister("ADC5", (adcMode == SDR)?0x403:0x401);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Calibrates the AsAd ADC.
 * See email from N. Usher on July 13th, 2013.
 * See Texas Instruments ADS6422 Data Sheet.
 * See email from J. Pibernat on October 16th, 2013.
 */
void CoBoNode::adcCal(AdcMode adcMode)
{
	// Straight binary data format, DESKEW test pattern (serial stream of 1010..)
	asad()->writeRegister("ADC2", 0x2C0);
	// Custom test pattern = 0xC3 ?
	asad()->writeRegister("ADC3", 0x0C3);
	asad()->writeRegister("ADC4", 0x000);
	// 2 Wire interface, bit clock selection
	asad()->writeRegister("ADC5", (adcMode == SDR)?0x403:0x401);

	// Initiate the CoBo deserializer deskew sequence
	ctrl()->writeField("pipeCtrl", "adcDeskew", 1);
	sleep_ms(1000);
	ctrl()->writeField("pipeCtrl", "adcDeskew", 0);

	// Straight binary data format, CUSTOM test pattern (registers ADC3 and ADC4)
	asad()->writeRegister("ADC2", 0x2A0);

	// Initiate the CoBo deserializer calibration sequence
	ctrl()->writeField("pipeCtrl", "bitSlipCalibrate", 1);
	sleep_ms(1000);
	ctrl()->writeField("pipeCtrl", "bitSlipCalibrate", 0);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configures the activation of the AsAd ADC channels and of the readout and sampling commands.
 * @see Texas Instruments ADS6422 Data Sheet.
 * @see AsAd Slow Control Protocol & Registers Mapping
 * @see AsAd Switched Capacitor Array Management, "AGET Inputs Manager".
 *
 * - ADC register A (ADC0):
 * 	 - D10	\<RST\> 		S/W RESET
 * 	 - D9	0
 * 	 - D8	0
 * 	 - D7	0
 * 	 - D6	0
 * 	 - D5	\<REF\> 		INTERNAL OR EXTERNAL
 * 	 - D4	\<PDN CHD\>	POWER DOWN CH D
 * 	 - D3	\<PDN CHC\>	POWER DOWN CH C
 * 	 - D2	\<PDN CHB\>	POWER DOWN CH B
 * 	 - D1	\<PDN CHA\>	POWER DOWN CH A
 * 	 - D0	\<PDN\>		GLOBAL POWER DOWN
 *
 * - Input Manager IM0 register: SWZ(3..0) Dbi(3..0)
 *   - SWZ(x) = 1	Connect ground to ZAP
 *   - SWZ(x) = 0	Connect VDD to ZAP
 *   - Dbi(x) = 1	Enable anti-sparkling protection for AGET
 *   - Dbi(x) = 0	Disable anti-sparkling protection for AGET
 *
 * - Input Manager IM1 register: ENW(3..0) ENR(3..0)
 *   - ENW(x) = 1	Enable WCK clock at the input of AGET x
 *                  Enable SCA Write signal to AGETx
 *   - ENW(x) = 0	Disable WCK clock at the input of AGET x
 *                  Disable SCA Write signal to AGETx e.g. disable sampling command
 *   - ENR(x) = 1	Enable RCK clock at the input of AGET x
 *                  Enable SCA Read signal to AGETx
 *   - ENR(x) = 0	Disable RCK clock at the input of AGET x
 *           	    Disable SCA Read signal to AGETx e.g. disable readout command
 */
void CoBoNode::asadAgetInputManagerInit(CCfg::CConfig nodeCfg)
{
	NODE_DEBUG() << "Configuring AsAd ADC channels and AsAd Inputs Manager";

	// Loop over AsAd boards
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
	{
		if (not isAsadActive(asadIdx)) continue;
		NODE_DEBUG() << "Activating AsAd no. " << asadIdx << " ADC channels";

		// Default register values
		uint8_t adc0 = 0;	// Enable all four ADC channels
		uint8_t im0 = 0xF; 	// Connect VDD to ZAP, Enable anti-sparkling protection
		uint8_t im1 = 0xFF;	// Allow sampling and readout for all four AGET chips

		// Loop over AGET chips
		for (size_t agetIdx = 0; agetIdx < 4; ++agetIdx)
		{

			// Check whether to activate slot
			bool isChipActive = isAgetActive(asadIdx, agetIdx);
			bool enablePosition = true;
			bool connectGroundToZAP = false; // SWZ
			bool antiSparklingEnabled = true; // Dbi

			try
			{
				enablePosition = agetConfig(nodeCfg, asadIdx, agetIdx)("Control")("enableAsAdPosition").getBoolValue();
				connectGroundToZAP = agetConfig(nodeCfg, asadIdx, agetIdx)("Control")("connectGroundToZAP").getBoolValue();
				antiSparklingEnabled = agetConfig(nodeCfg, asadIdx, agetIdx)("Control")("antiSparklingEnabled").getBoolValue();
			}
			catch (const CCfg::CConfig::NotFound &)
			{}

			// Disable slot
			if (not isChipActive or not enablePosition)
			{
				NODE_DEBUG() << "Disabling position no. " << agetIdx << " of AsAd no. " << asadIdx;
				// Power down ADC channel
				::utl::BitFieldHelper< uint8_t >::setBit(adc0, 1 + agetIdx, true);
				// Disable sampling command on AGET input
				::utl::BitFieldHelper< uint8_t >::setBit(im1, 4 + agetIdx, false);
				// Disable readout command on AGET input
				::utl::BitFieldHelper< uint8_t >::setBit(im1, agetIdx, false);
			}

			// Connect/disconnect VDD to/from ZAP interface
			BitFieldHelper< uint8_t >::setBit(im0, 4 + agetIdx, connectGroundToZAP);
			// Enable/disable anti-sparkling protection for AGET
			BitFieldHelper< uint8_t >::setBit(im0, agetIdx, antiSparklingEnabled);
		}

		// Set mask for broadcast of AsAd slow control command
		selectAsAdForSlowControl(asadIdx);

		// Send command
		asadInputsManagerConfig(adc0, im0, im1);
	}

	setAsAdSlowControlMask(asadMask_);
}
//_________________________________________________________________________________________________
/**
 * Set AGET input protections.
 * Two goals:
 * 	- Enable AsAd ZAP protections since AsAs is now totally powered.
 * 	- Prevent AGET sampling or readout commands until clocks are prepared.
 *
 */
void CoBoNode::asadSetAgetInputProtections()
{
	NODE_DEBUG() << "Setting AGET input protections";

	// Write to all AsAd boards simultaneously
	setAsAdSlowControlMask(asadMask_);

	// Configure connection to TPC pads and anti-sparkling protection
	asad()->writeRegister("IM0", 0x0F);
	// Disable sampling and readout commands on some AGET inputs
	asad()->writeRegister("IM1", 0x00);
}
//_________________________________________________________________________________________________
void CoBoNode::asadPowerStatus(const size_t & asadIdx, bool & plg, bool & al, bool & pwu)
{
	plg = isAsAdConnected(asadIdx);
	al = isAsAdPowerSupplied(asadIdx);
	pwu = isAsAdPoweredOn(asadIdx);
}
//__________________________________________________________________________________________________
/**
 * Reads the four 8-bit registers of the AsAd Identity Manager device and returns their value.
 *
 * Identifier: \<AIM0\>\<AIM2\>\<AIM1\>\<AIM3\>
 * These registers need to be read in the right order.
 * The AIM device is the fastest one to read, while the AsAd readout delay calibration sequence is based on the slowest.
 *
 * The AsAd identifier is a four 8-bit register (or a 32-bit register, if you prefer) composed of
 * - three 8-bit "read -only" fields (serial number etc...)
 * - and one 8-bit "read-write" field that allows the user to define its own label to identify a particular AsAd.
 * @see GET-AS-002-010: AsAd Slow Control Manager: Protocol & Registers Mapping
 */
uint32_t CoBoNode::asadGetIdentifier()
{
	uint32_t identifier = 0;
	try
	{
		// Read AIM registers one by one, in the right order
		for (size_t i=0; i < 4; ++i)
		{
			const std::string regName = std::string("AIM") + boost::lexical_cast< std::string >(i);
			const uint8_t regValue = asad()->readRegister(regName);
			utl::BitFieldHelper< uint32_t >::setField(identifier, 8u*(3u-i), 8u, regValue);
		}
	}
	catch (mdaq::utl::CmdException & e)
	{
		;
	}
	return identifier;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets AsAd board user ID.
 *
 * The AsAd identifier is a four 8-bit register (or a 32-bit register, if you prefer) composed of
 * - three 8-bit "read -only" fields (serial number etc...)
 * - and one 8-bit "read-write" field that allows the user to define its own label to identify a particular AsAd.
 */
void CoBoNode::asadSetBoardUserID(const uint8_t & userID)
{
	asad()->writeRegister("AIM0", userID);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Applies given configuration for AGET positions in AsAd.
 * Sets first register of the AsAd ADC device and the two registers of the AsAd Inputs Manager.
 *
 * See AsAd Switched Capacitor Array Management, "AGET Inputs Manager".
 * See AsAd Slow Control Manager: Protocol & registers mapping.
 *
 * Input Manager 'IM0' register: SWZ (3 .. 0) Dbi(3 .. 0)
 */
void CoBoNode::asadInputsManagerConfig(const uint8_t & adc0, const uint8_t & IM0, const uint8_t & IM1)
{
	NODE_DEBUG() << "Setting ADC0 to " << std::showbase << std::hex << (uint16_t) adc0 << ", IM0 to "
			<< (uint16_t) IM0 << " and IM1 to " << (uint16_t) IM1 << std::dec;
	// Power down some AGET's ADC channels
	asad()->writeRegister("ADC0", adc0);
	// Configure connection to TPC pads and anti-sparkling protection
	asad()->writeRegister("IM0", IM0);
	// Disable sampling and readout commands on some AGET inputs
	asad()->writeRegister("IM1", IM1);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Reads a value from a register of the "aget" device.
 * @param asadId Index of the AsAd board.
 * @param agetId Index of the AGET within its AsAd container.
 * @param registerName Name of the register.
 *
 * It assumes the SC_En mask of the aget_sc_configuration register is already set.
 *
 * @return Returns the register value read.
 */
uint64_t CoBoNode::agetRegRSC(size_t asadId, size_t agetId, const std::string & registerName)
{
	// Select AGET to read from
	ctrl()->writeField("agetChipSelect", "d", 4 * asadId + agetId);
	// Read register value
	return aget()->readRegister(registerName);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Writes a value in a register of the "aget" device.
 * @param asadId Index of the AsAd board.
 * @param agetId Index of the AGET within its AsAd container.
 * @param registerName Name of the register.
 * @param value Value to write.
 * @param verbose Verbosity flag
 */
void CoBoNode::agetRegWSC(size_t asadId, size_t agetId, const std::string & registerName, const uint64_t & value, const bool verbose)
{
	if (verbose)
	{
		NODE_DEBUG() << "Writing value " << std::hex << value << std::dec << " to register " << registerName << " of AGET no. " << asadId*4 + agetId;
	}
	// Select AGET to write to
	selectAgetForSlowControl(asadId, agetId);
	// Write register value
	aget()->writeRegister(registerName, value);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Adjusts the 3-bit global threshold in AGET register 1.
 */
void CoBoNode::setGlobalRegThresholds(const uint32_t & value)
{
	NODE_DEBUG() << "Initializing AGET global DAC thresholds to " << value;

	for (size_t asadId = 0; asadId < 4; ++asadId)
	{
		for (size_t agetId = 0; agetId < 4; ++agetId)
		{
			if (isAgetActive(asadId, agetId))
			{
				selectAgetForSlowControl(asadId, agetId);
				setGlobalRegThreshold(asadId, agetId, value);
			}
		}
	}

	setAgetSlowControlMask(agetMask_);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Adjusts the 3-bit global threshold in AGET register 1.
 */
void CoBoNode::setGlobalRegThreshold(size_t asadId, size_t agetId, const uint32_t & value)
{
	//NODE_DEBUG() << "Setting global threshold for AGET no. " << asadId*4 + agetId << " to " << value;

	uint64_t regValue = agetRegRSC(asadId, agetId, "reg1");
	::utl::BitFieldHelper< uint64_t >::setField(regValue, 19, 3, value);

	agetRegWSC(asadId, agetId, "reg1", regValue, false);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configures automatic computation of SCR data delay.
 */
void CoBoNode::calibrateAutoReadDataDelay(CCfg::CConfig nodeCfg)
{
		bool calcEnabled = true;
		try
		{
			calcEnabled = nodeCfg("Module")("calibrateAutoReadDataDelay").getBoolValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			// Keep default value
		}
		if (calcEnabled)
			calibrateAutoReadDataDelay();
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Triggers the automatic calibration of the data delay
 *  (time from rising edge of SCR signal until data appears at the CoBo).
 *
 * The calibration must be done after the other cable delays are calibrated
 *  and the local clocks manager of the AsAd has been configured.
 * It should be done before any registers of the AGET are configured,
 *  since AGET reg1 and reg2 are modified as part of the calibration.
 *
 * You can choose between the automatic or manual value by setting the SCRAutoDelayEnable# fields in the pipeCtrl register.
 */
void CoBoNode::calibrateAutoReadDataDelay()
{
	NODE_DEBUG() << "Calibrating SCR data delays";
	// Based on script dataCal.ecc
	ctrl()->writeField("aget_sc_configuration", "SC_En", agetMask_);

	aget()->writeRegister("reg1", 0x3e3a41a7);
	aget()->writeRegister("reg2", 0x0028a197);

	try
	{
		ctrl()->writeField("pipeCtrl", "SCRAutoDelayEnable", 1);
	}
	catch (const CmdException &)
	{
		ctrl()->writeField("pipeCtrl", "SCRAutoDelayEnable0", 1); // Field renamed on April 11th, 2013
	}

	ctrl()->writeField("pipeCtrl", "SCRAutoDelayCalibrate", 1);
	pipeInterruptEnable(false);
	ctrl()->writeField("pipeCtrl", "SCR", 0);
	ctrl()->writeField("pipeCtrl", "SCW", 0);
	ctrl()->writeField("pipeCtrl", "ScrScwOverride", 1);
	ctrl()->writeField("pipeCtrl", "reset", 0);
	ctrl()->writeField("pipeCtrl", "firmwareEnable", 1);
	ctrl()->writeField("pipeCtrl", "SCW", 1);
	sleep_ms(200);
	ctrl()->writeField("pipeCtrl", "SCW", 0);
	ctrl()->writeField("pipeCtrl", "SCR", 1);
	ctrl()->writeField("pipeCtrl", "SCR", 0);
	sleep_ms(200);
	ctrl()->writeField("pipeCtrl", "firmwareEnable", 0);
	ctrl()->writeField("pipeCtrl", "SCRAutoDelayCalibrate", 0);
	ctrl()->writeField("pipeCtrl", "ScrScwOverride", 0);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Selects automatic or manual SCR data delay, depending on the value of the
 * Setup[Conditions].CoBo[*].AsAd[*].Control.enableAutoReadDataDelay parameter or
 * Setup[Conditions].CoBo[*].Module.enableAutoReadDataDelay parameter
 * in the the configuration.
 * @param nodeCfg Configuration for the board.
 */
void CoBoNode::selectReadDataDelay(CCfg::CConfig nodeCfg)
{
	for (size_t asadId = 0; asadId < 4u; asadId++ )
	{
		asadSetAutoReadDataDelayEnabled(asadId, isSCRAutoDelayEnabled(nodeCfg, asadId));
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Enables or disables the use of automatically calibrated read data delays.
 *
 * Sets the SCRAutoDelayEnable CoBo register.
 * FIXME: Ignored in the microTCA CoBo firmware?
 * See email from Nathan Usher on April 11th and August 2nd 2013.
 */
void CoBoNode::asadSetAutoReadDataDelayEnabled(size_t asadId, bool enabled)
{
	if (isAsadActive(asadId))
	{
		NODE_DEBUG() << "Selecting " << (enabled?"automatic":"manual") << " data delays for AsAd no. " << asadId;
		std::string fieldName = std::string("SCRAutoDelayEnable");
		try
		{
			ctrl()->writeField("pipeCtrl", fieldName, enabled?1:0);
		}
		catch (const CmdException &)
		{
			fieldName = std::string("SCRAutoDelayEnable") + boost::lexical_cast<std::string>(asadId);
			ctrl()->writeField("pipeCtrl", "SCRAutoDelayEnable0", enabled?1:0); // Field renamed on April 11th, 2013
		}
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configure and start AsAd monitoring.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring
 * See email from J. Pibernat, Apr 7, 2017:
 * The goal of this step is to ensure the AsAd boards are fully powered
 * (the previous steps only showed that the FPGA and the monitoring chips were powered).
 * Since each GET setup uses a different power supply, the powering delay can vary a lot from setup to setup.
 * Therefore, we should ensure that the VDD voltage is high enough and, for that,
 * we need to first configure and start the monitoring device.
 * At this step, we only want to enable the monitoring, not the protection of the board.
 * So, before starting the monitoring we disable auto-protection.
 */
void CoBoNode::asadMonitoringInit(CCfg::CConfig nodeCfg)
{
	NODE_DEBUG() << "Initializing AsAd monitoring registers";
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
	{
		if (isAsadActive(asadIdx))
		{
			// Select AsAd board for slow control
			selectAsAdForSlowControl(asadIdx);

			// Reset AsAd monitoring device
			asadMonitoringReset();

			// Disable auto-protection
			asadDisableAutoProtection(asadIdx);

			// Define upper and lower limits
			asadMonitoringSetLimits(nodeCfg, asadIdx);

			// Define operational mode
			const bool fastMode = asadConfig(nodeCfg, asadIdx)("Monitoring")("fastConversionEnabled").getBoolValue(false);
			const uint8_t singleChannelItem = getAsAdMonitoringSingleChannel(nodeCfg, asadIdx);
			const bool singleChannelMode = asadConfig(nodeCfg, asadIdx)("Monitoring")("singleChannelModeEnabled").getBoolValue(false);
			const bool averagingDisabled = asadConfig(nodeCfg, asadIdx)("Monitoring")("averagingDisabled").getBoolValue(false);
			const bool enableAlarms = false;
			asadMonitoringSetOperationalMode(enableAlarms, fastMode, singleChannelItem, singleChannelMode, averagingDisabled);
		}
	}
	setAsAdSlowControlMask(asadMask_);
}
//__________________________________________________________________________________________________
/**
 * Initializes AsAd monitoring registers.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring
 */
void CoBoNode::asadAlertsInit(CCfg::CConfig nodeCfg)
{
	NODE_DEBUG() << "Configuring AsAd auto-protection...";
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
	{
		if (isAsadActive(asadIdx))
		{
			// Check AsAd board is not already in alarm
			if (asadAlarmBit(asadIdx) != 0)
			{
				std::ostringstream msg;
				msg << "AsAd board no. " << asadIdx << " is already in alarm mode!";
				throw mdaq::utl::CmdException(msg.str());
			}

			// Select AsAd board for slow control
			selectAsAdForSlowControl(asadIdx);

			// Step 1: define upper and lower limits
			asadMonitoringSetLimits(nodeCfg, asadIdx);

			// Step 2: define auto-protection level
			asadSetAutoProtectionLevel(nodeCfg, asadIdx);

			// Step 3: define operational mode
			asadMonitoringSetOperationalMode(nodeCfg, asadIdx);
		}
	}
	setAsAdSlowControlMask(asadMask_);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Read AsAd voltages, currents, temperatures.
 */
void CoBoNode::asadMonitoringRead()
{
	setAsAdSlowControlMask(asadMask_);
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
	{
		if (isAsadActive(asadIdx))
		{
			// Read AsAd voltages, currents, temperatures
			readAsAdVdd(asadIdx);
			readAsAdVdd1(asadIdx);
			readAsAdIdd(asadIdx);
			const int16_t Ti = readAsAdIntTemp(asadIdx);
			NODE_DEBUG() << "Read AsAd " << asadIdx << " temperature far from ADC: " << Ti << " deg C";
			const int16_t Te = readAsAdExtTemp(asadIdx);
			NODE_DEBUG() << "Read AsAd " << asadIdx << " temperature close to ADC: " << Te << " deg C";
		}
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Resets AsAd monitoring device.
 */
void CoBoNode::asadMonitoringReset()
{
	//NODE_DEBUG() << "Resetting AsAd monitoring device";
	// Reset all registers to default value
	asad()->writeRegister("configuration2", 0x80);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Defines AsAd monitoring limits.
 */
void CoBoNode::asadMonitoringSetLimits(CCfg::CConfig & nodeCfg, size_t asadId)
{
	NODE_DEBUG() << "Defining AsAd no. " << asadId << " monitoring limits";
	// First disable AsAd monitoring and auto-protection to "avoid instabilities in the system"
	asad()->writeRegister("configuration1", 0x2c);

	// Set AsAd monitoring limits
	setAsAdExtTempLimits(nodeCfg, asadId);
	setAsAdVddLimits(nodeCfg, asadId);
	setAsAdIddLimits(nodeCfg, asadId);
	//setAsAdIntTempLimits(nodeCfg, asadId);
	//setAsAdVdd1Limits(nodeCfg, asadId);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Disable all criteria for entering auto-protection.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring.
 */
void CoBoNode::asadDisableAutoProtection(const size_t & asadIdx)
{
	const uint8_t interruptMask1 = 0xFF;
	NODE_DEBUG() << "Setting AsAd no. " << asadIdx << " auto-protection level to " << std::hex << std::showbase << (uint16_t) interruptMask1 << std::dec;
	asad()->writeRegister("interruptMask1", interruptMask1);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets AsAd board auto-protection level.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring
 */
void CoBoNode::asadSetAutoProtectionLevel(CCfg::CConfig & nodeCfg, size_t asadIdx)
{
	const bool enableAlarms = asadConfig(nodeCfg, asadIdx)("Monitoring")("isAlarmSet").getBoolValue();

	// Set either highest or lowest auto-protection level
	// "interruptMask1" (highest protection level: 0x23, lowest: 0xFF)
	//     B0: (reserved) Write 1 only (disables Tint high limit)
	//     B1: (reserved) Write 1 only (disables Tint low limit)
	//     B2: If set, disables auto-protection caused by Text high limit.
	//     B3: If set, disables auto-protection caused by Text low limit.
	//     B4: If set, disables auto-protection caused by Text out of order.
	//     B5: (reserved) Write 1 only.
	//     B6: If set, disables auto-protection caused by VDD limits.
	//     B7: If set, disables auto-protection caused by IDD limits.
	uint8_t interruptMask1 = enableAlarms?0x23:0xFF;
	if (enableAlarms)
	{
		bool flag = false;
		try
		{
			flag = asadConfig(nodeCfg, asadIdx)("Monitoring")("LimitsText")("disableAlertsTooHigh").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 2u, flag);
			flag = asadConfig(nodeCfg, asadIdx)("Monitoring")("LimitsText")("disableAlertsTooLow").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 3u, flag);
			flag = asadConfig(nodeCfg, asadIdx)("Monitoring")("LimitsText")("disableAlertsOutOfOrder").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 4u, flag);
			flag = asadConfig(nodeCfg, asadIdx)("Monitoring")("LimitsVDD")("disableAlerts").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 6u, flag);
			flag = asadConfig(nodeCfg, asadIdx)("Monitoring")("LimitsIDD")("disableAlerts").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 7u, flag);
		}
		catch (const CCfg::CConfig::NotFound & e)
		{
			NODE_WARN() << e.what();
		}
	}
	NODE_DEBUG() << "Setting AsAd no. " << asadIdx << " auto-protection level to " << std::hex << std::showbase << (uint16_t) interruptMask1 << std::dec;
	asad()->writeRegister("interruptMask1", interruptMask1);

	// Leaving "interruptMask2" to its default value (has to do with VDD1)
	// const uint8_t interruptMask2 = 0;
	// asad()->writeRegister("interruptMask2", interruptMask2);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets AsAd monitoring device operational mode.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring
 */
void CoBoNode::asadMonitoringSetOperationalMode(bool enableAlarms, bool averagingDisabled, bool singleChannelMode, uint8_t singleChannelItem, bool fastMode)
{
	// asad()->writeRegister("intTempOffset", 0x0);
	// asad()->writeRegister("extTempOffset", 0x0);

	// "configuration3" (default value:0x0)
	//                B0 0: ADC for monitoring in slow conversion mode (clock @  1.4 kHz)
	//                   1: ADC for monitoring in fast conversion mode (clock @ 22.5 kHz) and analog filters disabled
	//                B1 0: (not used)
	//                B2 0: (reserved)
	//                B3 0: (not used)
	//                B4 0: internal ADC reference given by an internal Vref
	//                B5 0: (not used)
	//                B6 0: (not used)
	//                B7 0: (reserved)
	uint8_t configuration3 = 0;
	::utl::BitFieldHelper< uint8_t>::setBit(configuration3, 0, fastMode);

	//NODE_DEBUG() << "Setting AsAd register 'configuration3' to " << std::hex << std::showbase << (uint16_t) configuration3 << std::dec;
	asad()->writeRegister("configuration3", configuration3);

	// "configuration2" (default value:0x0)
	//               <B0:B2> Defines the monitored item in single channel mode
	//                        000: VDD1; 001: Internal temperature; 010: External temperature; 100: VDD; 101: IDD
	//                        (Round robin mode has to be selected by default instead of single channel mode)
	//                B3 0: (reserved)
	//                B4 0: round robin mode
	//                   1: single channel mode
	//                B5 0: enable averaging every measurement on all item 16 times
	//                   1: disable averaging
	//                B6 0: (not used)
	//                B7 1: software reset
	uint8_t configuration2 = 0;
	::utl::BitFieldHelper< uint8_t>::setField(configuration2, 0, 3, singleChannelItem);
	::utl::BitFieldHelper< uint8_t>::setBit(configuration2, 4, singleChannelMode);
	::utl::BitFieldHelper< uint8_t>::setBit(configuration2, 5, averagingDisabled);

	//NODE_DEBUG() << "Setting AsAd register 'configuration2' to " << std::hex << std::showbase << (uint16_t) configuration2 << std::dec;
	asad()->writeRegister("configuration2", configuration2);

	// "configuration1" (default value:0x0)
	//    B0 0: stop monitoring
	//       1:start monitoring
	//    <B1:B2> 10: ADS6422 temperature measurement
	//    B3 1: VDD value measurement
	//    B4 0: (reserved)
	//    B5 0: Alert enabled
	//       1: Alert disabled
	//    B6 0: Alert output active low on AsAd (high on CoBo)
	//    PD 0: Power up
	//       1: Power down
	uint8_t configuration1 = 0x0D;
	::utl::BitFieldHelper< uint8_t>::setBit(configuration1, 5, enableAlarms?0:1);

	//NODE_DEBUG() << "Setting AsAd register 'configuration1' to " << std::hex << std::showbase << (uint16_t) configuration1 << std::dec;
	asad()->writeRegister("configuration1", configuration1);

	NODE_DEBUG() << "Set operational mode of AsAd monitoring device to:" << std::hex << std::showbase
			<< " cfg1=" << (uint16_t) configuration1
			<< " cfg2=" << (uint16_t) configuration2
			<< " cfg3=" << (uint16_t) configuration3 << std::dec;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets AsAd monitoring device operational mode.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring
 */
void CoBoNode::asadMonitoringSetOperationalMode(CCfg::CConfig & nodeCfg, size_t asadIdx)
{
	// "configuration3" (default value:0x0)
	//                B0 0: ADC for monitoring in slow conversion mode (clock @  1.4 kHz)
	//                   1: ADC for monitoring in fast conversion mode (clock @ 22.5 kHz) and analog filters disabled
	//                B1 0: (not used)
	//                B2 0: (reserved)
	//                B3 0: (not used)
	//                B4 0: internal ADC reference given by an internal Vref
	//                B5 0: (not used)
	//                B6 0: (not used)
	//                B7 0: (reserved)
	bool fastMode = asadConfig(nodeCfg, asadIdx)("Monitoring")("fastConversionEnabled").getBoolValue(false);

	// "configuration2" (default value:0x0)
	//               <B0:B2> Defines the monitored item in single channel mode
	//                        000: VDD1; 001: Internal temperature; 010: External temperature; 100: VDD; 101: IDD
	//                        (Round robin mode has to be selected by default instead of single channel mode)
	//                B3 0: (reserved)
	//                B4 0: round robin mode
	//                   1: single channel mode
	//                B5 0: enable averaging every measurement on all item 16 times
	//                   1: disable averaging
	//                B6 0: (not used)
	//                B7 1: software reset
	uint8_t singleChannelItem = getAsAdMonitoringSingleChannel(nodeCfg, asadIdx);
	bool singleChannelMode = asadConfig(nodeCfg, asadIdx)("Monitoring")("singleChannelModeEnabled").getBoolValue(false);
	bool averagingDisabled = asadConfig(nodeCfg, asadIdx)("Monitoring")("averagingDisabled").getBoolValue(false);

	// "configuration1" (default value:0x0)
	//    B0 0: stop monitoring
	//       1:start monitoring
	//    <B1:B2> 10: ADS6422 temperature measurement
	//    B3 1: VDD value measurement
	//    B4 0: (reserved)
	//    B5 0: Alert enabled
	//       1: Alert disabled
	//    B6 0: Alert output active low on AsAd (high on CoBo)
	//    PD 0: Power up
	//       1: Power down
	const bool enableAlarms = asadConfig(nodeCfg, asadIdx)("Monitoring")("isAlarmSet").getBoolValue();
	asadMonitoringSetOperationalMode(enableAlarms, fastMode, singleChannelItem, singleChannelMode, averagingDisabled);
}
// --------------------------------------------------------------------------------------------------------------------
void CoBoNode::asadMonitoringCheckVDD(CCfg::CConfig nodeCfg, const size_t maxReads, const double delay_ms)
{
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
	{
		if (isAsadActive(asadIdx))
		{
			bool ok = false;
			const double vddMin = asadConfig(nodeCfg, asadIdx)("Monitoring")("LimitsVDD")("lowerValue").getRealValue(3.2);
			NODE_INFO() << "Waiting for AsAd voltage in AsAd no. " << asadIdx << " to settle above " << vddMin << " V ...";
			for (size_t i=0; i < maxReads; ++i)
			{
				// Read VDD until high enough
				const double vdd = readAsAdVdd(asadIdx);
				if (vdd >= vddMin)
				{
					ok = true;
					break;
				}
				sleep_ms(delay_ms);
			}
			if (not ok)
			{
				std::ostringstream msg;
				msg << "VDD voltage of AsAd board no. " << asadIdx << " did not reach " << vddMin
						<< " V! Check the hardware setup of your power supply or contact CENBG for further information.";
				throw CmdException(msg.str());
			}

		}
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Enables or disables the systematic readout of Fixed Pattern Noise channels.
 */
void CoBoNode::coboEnableFPNRead(CCfg::CConfig & nodeCfg)
{
	// Compute flag value
	bool isFPNRead = false;
	// Loop over each of the 4 AsAd boards
	for (size_t asadId=0; asadId < 4u; ++asadId)
	{
		// Loop over AGET chips
		for (size_t agetId = 0; agetId < 4u; ++agetId)
		{
			if (not isAgetActive(asadId, agetId))
				continue;
			const bool isAgetFpnRead = agetConfig(nodeCfg, asadId, agetId)("Global")("Reg1")("isFPNRead").getBoolValue();
			isFPNRead = isFPNRead or isAgetFpnRead;
		}
	}

	// Set register value
	coboEnableFPNRead(isFPNRead);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Enables or disables the systematic readout of Fixed Pattern Noise channels.
 * Sets the corresponding CoBo register. The AGET register 1 also needs to be set accordingly.
 */
void CoBoNode::coboEnableFPNRead(bool isFPNRead)
{
	NODE_DEBUG() << "Setting flag for readout of FPN channels to " << isFPNRead;
	ctrl()->writeField("pipeCtrl", "isFPNRead", isFPNRead?1:0); // Always read four FPN channels, even in partial readout mode.
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets an AGET register to a given value, with an optional check.
 * @param asadId AsAd board identifier.
 * @param agetId AGET chip identifier within the AsAd board.
 * @param registerName Name of the AGET register to configure.
 * @param value Register value.
 * @param check Whether to read the register after setting its value, to check if the operation has been successful.
 */
void CoBoNode::agetSetRegister(size_t asadId, size_t agetId, const std::string & registerName, const uint64_t & value, bool check)
{
	// Write value
	agetRegWSC(asadId, agetId, registerName, value);

	// Check
	if (check)
	{
		// Read back value
		const uint64_t checkValue = agetRegRSC(asadId, agetId, registerName);
		if (value != checkValue)
		{
			std::ostringstream msg;
			msg << "Check of register '" << registerName << " of AGET " << agetId << "' of AsAd " << asadId << "' failed ("
					<< std::hex << checkValue << " instead of " << value << ").";
			// Handle special case of bit 21 of register 2 failure
			if ("reg2" == registerName)
			{
				uint64_t regXor =  value ^ checkValue;
				if (regXor == UINT64_C(0x200000))
				{
					NODE_WARN() << msg.str();
					return;
				}
			}
			throw CmdException(msg.str());
		}
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the AGET registers to the values deduced from the configuration.
 * @param nodeCfg Configuration to extract the register values from.
 * @param asadId AsAd board identifier.
 * @param agetId AGET chip identifier within the AsAd board.
 * @param check Whether to read the registers after setting their values to check if the operation has been successful.
 */
void CoBoNode::agetSetRegisters(CCfg::CConfig & nodeCfg, size_t asadId, size_t agetId, bool check)
{
	NODE_DEBUG() << "Configuring registers of AGET no. " << (4*asadId + agetId);
	if (isAgetActive(asadId, agetId))
	{
		// Get AGET chip version to customize configuration between AGET and ASTRE
		uint16_t chipVersion = agetVersion(asadId, agetId);

		const std::string regNames[13] = { "reg1", "reg2", "reg3", "reg4", "reg6", "reg7", "reg8", "reg8msb", "reg9", "reg9msb", "reg10", "reg11", "reg12" };
		for (size_t i=0; i < 13; ++i)
		{
			const uint64_t cfgValue = CoBoConfigParser::agetRegValue(nodeCfg, asadId, agetId, regNames[i], chipVersion);
			agetSetRegister(asadId, agetId, regNames[i], cfgValue, check);

			// For AGET reg12 (readout pointer offset), also set the corresponding CoBo register 'readDepth' (number of time buckets to read).
			if ("reg12" == regNames[i])
			{
				const uint16_t readoutPointerOffset = cfgValue;
				const bool is2pModeEnabled = nodeCfg("Module")("enableMem2pMode").getBoolValue();
				const uint16_t depth = is2pModeEnabled ? 256 - readoutPointerOffset : 512 - readoutPointerOffset;
				NODE_DEBUG() << "Setting readout depth to " << depth << " bucket(s) starting from offset " << readoutPointerOffset;
				ctrl()->writeRegister("readDepth", depth);
			}
		}
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets ReadIfHit registers for a given AGET chip according to values in configuration file.
 * @param nodeCfg Configuration.
 * @param asadIdx Index of the AsAd board the AGET belongs to (0 to 3).
 * @param agetIdx Index of the AGET chip within its AsAd board (0 to 3).
 */
void CoBoNode::agetSetReadIfHitMask(CCfg::CConfig & nodeCfg, size_t asadIdx, size_t agetIdx)
{
	if (isAgetActive(asadIdx, agetIdx))
	{
		// Select AGET chip
		ctrl()->writeRegister("maskSelect", asadIdx*4+agetIdx);
		// Get desired values
		const uint32_t val1_4 = agetReadoutMask(nodeCfg, asadIdx, agetIdx, "only_if_hit", 0, 3);
		const uint32_t val5_36 = agetReadoutMask(nodeCfg, asadIdx, agetIdx, "only_if_hit", 4, 35);
		const uint32_t val37_68 = agetReadoutMask(nodeCfg, asadIdx, agetIdx, "only_if_hit", 36, 67);
		// Write register values
		NODE_DEBUG() << "Setting 'readIfHitMask' registers for AGET no. " << asadIdx*4+agetIdx << " to " << std::showbase << std::hex << val1_4 << ", " << val5_36 << ", " << val37_68 << std::dec;
		ctrl()->writeRegister("readIfHitMask1_4", val1_4);
		ctrl()->writeRegister("readIfHitMask5_36", val5_36);
		ctrl()->writeRegister("readIfHitMask37_68", val37_68);
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets ReadAlways registers for a given AGET chip according to values in configuration file.
 * @param nodeCfg Configuration.
 * @param asadIdx Index of the AsAd board the AGET belongs to (0 to 3).
 * @param agetIdx Index of the AGET chip within its AsAd board (0 to 3).
 */
void CoBoNode::agetSetReadAlwaysMask(CCfg::CConfig & nodeCfg, size_t asadIdx, size_t agetIdx)
{
	if (isAgetActive(asadIdx, agetIdx))
	{

		// Get desired values
		const uint32_t val1_4 = agetReadoutMask(nodeCfg, asadIdx, agetIdx, "always", 0, 3);
		const uint32_t val5_36 = agetReadoutMask(nodeCfg, asadIdx, agetIdx, "always", 4, 35);
		const uint32_t val37_68 = agetReadoutMask(nodeCfg, asadIdx, agetIdx, "always", 36, 67);

		NODE_DEBUG() << "Setting 'readAlwaysMask' registers for AGET no. "
			<< asadIdx*4+agetIdx << " to " << std::showbase << std::hex
			<< val1_4 << ", "
			<< val5_36 << ", " <<
			val37_68 << std::dec;

		hwNode().regProc("ctrl")
			.write("maskSelect", asadIdx*4+agetIdx) // Select AGET chip
			.write("readAlwaysMask1_4", val1_4)
			.write("readAlwaysMask5_36", val5_36)
			.write("readAlwaysMask37_68", val37_68);
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Switches to fast control.
 * This function switches the state machine of the CoBo AGET controller
 * and the state machine on the AGET to the fast control mode,
 * i.e. the reading/writing of the hit register
 */
void CoBoNode::switchToFastControl(uint16_t agetMask)
{
	hwNode().regProc("ctrl")
		// Resets the state machine for the AGET slow control.
		.write("agetController", "reset", 1) // Will reset the error flags and the busy flag (needs to be pulsed: 0,1,0).
		.write("agetController", "reset", 0)
		// Makes a fake slow control to send the three pulses
		// necessary for AGET to switch to the fast control mode
		.write("aget_sc_configuration", "SC_En", 0) // Mask indicating which of the 16 AGet's the slow control commands should be broadcast to.
		.write("aget_sc_configuration", "pulseEnable", agetMask) // Needed when switching between slow control mode and hit register mode. Bits written into the first Byte read by the FPGA.
		.write("agetController", "hitRegister", 1) // If set, switches to hit register mode.
		.write("agetController", "calibration", 0) // If set, calibration delays are automatically computed.
		.write("aget_sc_rw", 1) // Triggers the read (1) or write (0) operation of the AGET register.
		.sleep_ms(100)
		// Enables the control of the AGET and disables the pulse generation
		// No 3-pulse should be sent during a fast control of the AGET
		.write("aget_sc_configuration", "SC_En", agetMask)
		.write("aget_sc_configuration", "pulseEnable", 0);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Switches to fast control.
 * This function switches the state machine of the CoBo AGET controller
 * and the state machine on the AGET to the fast control mode,
 * i.e. the reading/writing of the hit register
 */
void CoBoNode::switchToFastControl()
{
	switchToFastControl(agetMask_);
}
// --------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------
/**
 * Switches to slow control.
 * This function switches the state machine of the CoBo AGET controller
 * and the state machine of the AGET to the slow control mode
 */
void CoBoNode::switchToSlowControl(uint16_t agetMask)
{
	hwNode().regProc("ctrl")
		// Resets state machine for AGET slow control and disables firmware hand on the AGET slow control
		.write("pipeCtrl", "firmwareEnable", 0)
		.write("agetController", "reset", 1)
		.write("agetController", "reset", 0)
		// Issues fake slow control to send the 3 pulses necessary for AGET to switch to the slow control mode
		.write("aget_sc_configuration", "SC_En", 0)
		.write("aget_sc_configuration", "pulseEnable", agetMask)
		.write("agetController", "hitRegister", 0)
		.write("agetController", "calibration", 0)
		.write("aget_sc_rw", 1)
		.sleep_ms(1000)
		// Enables AGET control and disables pulse generation: no 3-pulse should be sent during fast control
		.write("aget_sc_configuration", "SC_En", agetMask)
		.write("aget_sc_configuration", "pulseEnable", 0);
}

// --------------------------------------------------------------------------------------------------------------------
/**
 * Switches to slow control.
 * This function switches the state machine of the CoBo AGET controller
 * and the state machine of the AGET to the slow control mode
 */
void CoBoNode::switchToSlowControl()
{
	switchToSlowControl(agetMask_);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Selects, on each specified AGET chip, a single channel to always read.
 * It is assumed that we already are in slow control mode.
 */
void CoBoNode::selectSingleAgetReadoutChannel(uint16_t agetMask, size_t chanId)
{
	NODE_DEBUG() << "Selecting single channel (no. " << chanId << ") per chip to always read";

	RemoteRegisterProcessor& regProc = hwNode().regProc("ctrl");

	// Loop over AGET chips of all AsAd boards
	for (size_t agetId=0; agetId < 16; ++agetId)
	{
		// Skip inactive chips
		if (not utl::BitFieldHelper< uint16_t >::getField(agetMask, agetId, 1)) continue;

		regProc	<< Cmd<Write>("maskSelect", "", agetId) // Tell CoBo which chip we are setting mask for
				<< Cmd<Write>("readIfHitMask1_4", "", 0) // Set mask of channels to read only if hit
				<< Cmd<Write>("readIfHitMask5_36", "", 0)
				<< Cmd<Write>("readIfHitMask37_68", "", 0);

		// Set mask of channels to always read
		if (chanId < 4)
		{
			regProc	<< Cmd<Write>("readAlwaysMask1_4", "", UINT32_C(1) << (3 - chanId))
					<< Cmd<Write>("readAlwaysMask5_36", "", 0)
					<< Cmd<Write>("readAlwaysMask37_68", "", 0);
		}
		else if (chanId < 36)
		{
			regProc	<< Cmd<Write>("readAlwaysMask1_4", "", 0)
					<< Cmd<Write>("readAlwaysMask5_36", "", UINT32_C(1) << (35 - chanId))
					<< Cmd<Write>("readAlwaysMask37_68", "", 0);
		}
		else if (chanId < 68)
		{
			regProc	<< Cmd<Write>("readAlwaysMask1_4", "", 0)
					<< Cmd<Write>("readAlwaysMask5_36", "", 0)
					<< Cmd<Write>("readAlwaysMask37_68", "", UINT32_C(1) << (67 - chanId));
		}
		else
		{
			regProc	<< Cmd<Write>("readAlwaysMask1_4", "", 0)
					<< Cmd<Write>("readAlwaysMask5_36", "", 0)
					<< Cmd<Write>("readAlwaysMask37_68", "", 0);
		}
	}
	regProc.flush();
}
void CoBoNode::selectSingleAgetReadoutChannel(size_t chanId)
{
	selectSingleAgetReadoutChannel(agetMask_, chanId);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Selects a single test channel per chip.
 * It is assumed that we already are in slow control mode.
 */
void CoBoNode::selectSingleAgetTestChannel(size_t chanId)
{
	NODE_DEBUG() << "Selecting single channel (no. " << chanId << ") per chip to test";

	RemoteRegisterProcessor& regProc = hwNode().regProc("aget");
	if (chanId < 34)
	{
		regProc	<< Cmd<Write>("reg3", "", UINT64_C(1) << (33 - chanId))
				<< Cmd<Write>("reg4", "", 0);
	}
	else if (chanId < 68)
	{
		regProc	<< Cmd<Write>("reg3", "", 0)
				<< Cmd<Write>("reg4", "", UINT64_C(1) << (chanId - 34));
	}
	else
	{
		regProc	<< Cmd<Write>("reg3", "", 0)
				<< Cmd<Write>("reg4", "", 0);
	}
	regProc.flush();
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets LSB threshold of an AGET channel.
 * It is assumed that we already are in slow control mode.
 */
void CoBoNode::setAgetLsbThreshold(const size_t & channelId, const uint8_t threshold)
{
	RemoteRegisterProcessor& regProc = hwNode().regProc("aget");
	if (channelId < 68)
	{
		NODE_DEBUG() << "Setting LSB threshold of AGET channel no. " << std::dec << channelId << " to " << uint16_t(threshold);
		// Find register
		size_t physicalChannelId = 1 + CoBoConfigParser::agetPhysicalChannelId(channelId);
		if (physicalChannelId <= 0)
			return;
		std::string regName = "reg8";
		if (physicalChannelId >= 49)
			regName = "reg9msb";
		else if (physicalChannelId >= 33)
			regName = "reg9";
		else if (physicalChannelId <= 16)
			regName = "reg8msb";
		// Read current register value (from one specific chip!)
		uint64_t regValue;
		regProc.read(regName, (Ice::Long &) regValue);
		// NODE_DEBUG() << "Old value: " << std::hex << std::showbase << regValue;
		// Compute new register value
		size_t offset = agetChannelLsbThresholdOffset(regName, channelId);
		// NODE_DEBUG() << "Offset = " << offset;
		::utl::BitFieldHelper<uint64_t>::setField(regValue, offset, 4, threshold);
		// Write new register value (to all chips!)
		// NODE_DEBUG() << "New value: " << std::hex << std::showbase << regValue;
		regProc.write(regName, regValue);
	}
	else
	{
		NODE_DEBUG() << "Setting LSB threshold of all AGET channels to " << uint16_t(threshold);
		// Build new register value
		uint64_t regValue = 0;
		for (size_t offset=0; offset < 64; offset += 4)
		{
			::utl::BitFieldHelper<uint64_t>::setField(regValue, offset, 4, threshold);
		}
		// NODE_DEBUG() << "Writing " << std::hex << std::showbase << regValue;
		// Write to all concerned registers
		const std::string regNames[4] = {"reg8", "reg8msb", "reg9", "reg9msb"};
		for (size_t ireg=0; ireg < 4; ++ireg)
		{
			regProc.write(regNames[ireg], regValue);
		}
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Enables or disables the modification of the hit pattern before readout.
 * If enabled, the hit pattern is OR'ed and AND'ed with registers 'readAlways' and 'readIfHit' respectively.
 * @param nodeCfg Configuration.
 */
void CoBoNode::setWriteHitRegisterMask(CCfg::CConfig nodeCfg)
{
	unsigned enableBit = nodeCfg("Module")("enableWriteHittedregister").getBoolValue() ? 1 : 0;
	NODE_DEBUG() << "Modification of hit pattern before readout has been " << (enableBit == 1 ? "enabled" : "disabled");
	ctrl()->writeField("pipeCtrl", "writeHitRegisterMaskEnable", enableBit);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Enables or disables the test mode, where frames are filled with an incremental counter instead of actual data.
 */
void CoBoNode::setTestMode(CCfg::CConfig & nodeCfg)
{
	// Get desired value
	unsigned enableTestModeBit = nodeCfg("Module")("isAcqTestMode").getBoolValue() ? 1 : 0;
	// Write
	ctrl()->writeField("pipeCtrl", "testMode", enableTestModeBit);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configures the subtraction of the moving average and deviation of each AGET's multiplicity signal.
 *
 * To compensate for the noise on the AGET multiplicity output, the CoBo now calculates a moving average
 *  and moving absolute deviation of the multiplicity signal,
 *  and those values can be subtracted from the multiplicity signal to shift the baseline down to 0.
 */
void CoBoNode::setCoBoMultSubtractMode(CCfg::CConfig & nodeCfg)
{
	// Get desired values from config
	bool subtractAverage = false;
	uint8_t subtractDevs = 0;
	try
	{
		subtractAverage = nodeCfg("Module")("multSubtractAverage").getBoolValue();
		subtractDevs = nodeCfg("Module")("multSubtractDevs").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		// keep default value
	}

	// Set register fields
	try
	{
		NODE_DEBUG() << "Setting multiplicity subtractAverage to " << subtractAverage << " and subtractDevs to " <<  (short) subtractDevs;

		hwNode().regProc("ctrl")
			// When set, average multiplicity value is subtracted from AGET output to shift baseline down to 0.
			.write("pipeCtrl", "MultSubtractAverage", subtractAverage ? 1 : 0)
			// Twice the number of absolute deviations to subtract from multiplicity signal.
			// NOTE: A value of 4 (2 deviations) worked well to prevent noise from causing triggers.
			.write("pipeCtrl", "MultSubtractDevs", subtractDevs);
	}
	catch(const CmdException& e)
	{
		NODE_WARN() << "Error setting external trigger mode: " << e.reason;
		return;
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configures the 4 bits of the CoBo 'pipeCtrl' register indicating which of the 4 chips to ignore in the multiplicity trigger mode.
 *
 * On Micro-TCA CoBo, the same setting applies to all AsAd boards.
 */
void CoBoNode::setCoBoSuppressMultiplicity(CCfg::CConfig & nodeCfg)
{
	// Loop over AGET positions
	uint8_t suppressMultiplicity = 0;
	for (size_t agetIdx = 0; agetIdx < 4; agetIdx++)
	{
		// Check in configuration file
		size_t suppressChipCount = 0;
		size_t activeChipCount = 0;
		for (size_t asadIdx = 0; asadIdx < 4u; ++asadIdx)
		{
			if (not isAsadActive(asadIdx)) continue;


			if (isAgetActive(asadIdx, agetIdx))
			{
				try
				{
					activeChipCount++;
					if (agetConfig(nodeCfg, asadIdx, agetIdx)("Global")("suppressMultiplicity").getBoolValue())
					{
						suppressChipCount++;
					}

				}
				catch (const CCfg::CConfig::NotFound &)
				{
					// keep default value
				}
			}
		}
		// Check consistency
		bool suppressChipPosition = false;
		if (activeChipCount <= 0)
		{
			suppressChipPosition = true;
		}
		else if (suppressChipCount > 0)
		{
			if (suppressChipCount == activeChipCount)
				suppressChipPosition = true;
			else if (suppressChipCount < activeChipCount)
			{
				std::ostringstream msg;
				msg << "Cannot suppress multiplicity from AGET position " << agetIdx << " on one of CoBo '" << id().toString()
								<< "' 's AsAd board while keeping it on other AsAd boards. Check parameter 'suppressMultiplicity'";
				throw ::mdaq::utl::CmdException(msg.str());
			}
		}
		// Suppress multiplicity
		utl::BitFieldHelper<uint8_t>::setField(suppressMultiplicity, agetIdx, 1, suppressChipPosition ? 1 : 0);
	}

	// Set CoBo register
	NODE_DEBUG() << "Setting 'suppressMultiplicity' field to " << std::hex << std::showbase << (unsigned short) suppressMultiplicity << std::noshowbase << std::dec;
	ctrl()->writeField("pipeCtrl", "suppressMultiplicity", suppressMultiplicity);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Enables zero suppression and sets the polarity of the zero suppression threshold.
 * @param nodeCfg Configuration
 */
void CoBoNode::setCoBoZeroSuppressionMode(CCfg::CConfig& nodeCfg)
{
	// Enable or disable zero suppression
	// Get desired value from configuration
	bool enableSuppression = nodeCfg("Module")("enableZeroSuppression").getBoolValue(false);
	NODE_DEBUG() << (enableSuppression ? "Enabling" : "Disabling") << " zero suppression";


	// Set the polarity of the zero suppression threshold.
	bool zeroSuppressInvert = nodeCfg("Module")("zeroSuppressionInverted").getBoolValue(false);
	if (enableSuppression and zeroSuppressInvert)
	{
		NODE_DEBUG() << "CoBo will suppress sample values ABOVE zero suppression thresholds.";
	}
	try
	{
		hwNode().regProc("ctrl")
			.write("pipeCtrl", "zeroSuppressionEnable", enableSuppression ? 1 : 0)
			.write("pipeCtrl", "zeroSuppressInvert", zeroSuppressInvert ? 1 : 0);
	}
	catch (const CmdException & e)
	{
		NODE_WARN() << e.what(); // This bit was described after email from N. U. on May 20th, 2014.
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets CoBo registers for zero suppression thresholds for a given AGET chip according to values in configuration file.
 * @param nodeCfg Configuration.
 * @param asadIdx Index of the AsAd board the AGET belongs to (0 to 3).
 * @param agetIdx Index of the AGET chip within its AsAd board (0 to 3).
 */
void CoBoNode::setCoBoZeroSuppressionThresholds(CCfg::CConfig & nodeCfg, size_t asadIdx, size_t agetIdx)
{
	if (isAgetActive(asadIdx, agetIdx))
	{
		NODE_DEBUG() << "Setting zero suppression thresholds for AGET no. " << asadIdx*4+agetIdx;
		try
		{
			for (size_t channelId=0; channelId < 68; ++channelId)
			{
				// Get threshold value from configuration
				const uint16_t threshold = agetConfig(nodeCfg, asadIdx, agetIdx)("channel", channelId)("zeroSuppressionThreshold").getIntValue();

				// Set register
				std::ostringstream regName;
				regName << "aget" << agetIdx << "ch" << channelId;
				zeroSuppressDevice(asadIdx)->writeField(regName.str(), "threshold", threshold);
			}
		}
		catch (const CCfg::CConfig::NotFound & e)
		{
			NODE_WARN() << e.what() << " Zero suppression thresholds have not been configured.";
		}
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Enables or disables the 2p mode, where the SCA is split in two separate memories of 256 cells.
 */
void CoBoNode::enableMem2pMode(CCfg::CConfig & nodeCfg)
{
	bool enable2pMode = nodeCfg("Module")("enableMem2pMode").getBoolValue();
	NODE_DEBUG() << "Setting 'modeEnable_2p' to " << enable2pMode;
	ctrl()->writeField("pipeCtrl", "modeEnable_2p", enable2pMode?1:0);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configures the value of the 'dataSource' field in the header of the data frames.
 * This Data source ID typically refers to the detector or detector type from which the frame data have been extracted.
 * FIXME: We should ensure that all boards use the same dataSource value.
 */
void CoBoNode::setDataSourceId(CCfg::CConfig & nodeCfg)
{
	uint8_t dataSourceId = 0;
	try
	{
		dataSourceId = nodeCfg("Module")("dataSource").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		NODE_WARN() << e.what();
	}
	ctrl()->writeField("dataSource", "value", dataSourceId);
	NODE_DEBUG() << "Setting frame data source ID to " << (short) dataSourceId;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configures the value of the CoBo index in the header of the data frames.
 * FIXME: We should ensure that each board uses a different ID value.
 */
void CoBoNode::setCoBoId(CCfg::CConfig & nodeCfg)
{
	uint8_t coboIdentifier = 0;
	try
	{
		coboIdentifier = nodeCfg("Module")("coboId").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		NODE_WARN() << e.what();
	}
	ctrl()->writeField("COBO_index", "d", coboIdentifier);
	NODE_DEBUG() << "Setting frame CoBo ID to " << (short) coboIdentifier;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets period of periodic trigger.
 */
void CoBoNode::setTriggerPeriod(CCfg::CConfig & nodeCfg)
{
	uint32_t period_10ns = 100000000u; // 1s
	try
	{
		period_10ns = nodeCfg("Module")("triggerPeriod_10ns").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		period_10ns = 100000*nodeCfg("Module")("triggerPeriod").getIntValue(); // ms
	}
	ctrl()->writeRegister("triggerPeriod_10ns", period_10ns);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets delay between event causing the trigger and the actual freeze of the ASIC memory.
 */
void CoBoNode::setTriggerDelay(CCfg::CConfig & nodeCfg)
{
	int delay_10ns = nodeCfg("Module")("triggerDelay").getIntValue();
	ctrl()->writeRegister("triggerDelay", delay_10ns);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets delay between event causing the second trigger and the actual freeze of the ASIC memory.
 */
void CoBoNode::setTriggerDelay2p(CCfg::CConfig & nodeCfg)
{
	int delay_10ns = nodeCfg("Module")("triggerDelay_2p").getIntValue();
	ctrl()->writeRegister("triggerDelay_2p", delay_10ns);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the delay after the first trigger after which the second trigger can be issued.
 */
void CoBoNode::setTriggerDeadTime2p(CCfg::CConfig & nodeCfg)
{
	int deadTime2p_10ns = nodeCfg("Module")("triggerDeadTime_2p").getIntValue();
	ctrl()->writeRegister("triggerWindowStart_2p", deadTime2p_10ns);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the trigger period of the second trigger (if periodic).
 */
void CoBoNode::setTriggerTime2p(CCfg::CConfig & nodeCfg)
{
	int time2p_10ns = nodeCfg("Module")("triggerTime_2p").getIntValue();
	ctrl()->writeRegister("triggerTime_10ns_2p", time2p_10ns);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the maximum duration to wait for the second trigger.
 *
 * After this delay, the first trigger is canceled.
 */
void CoBoNode::setTriggerTimeOut2p(CCfg::CConfig & nodeCfg)
{
	int timeOut2p_10ns = nodeCfg("Module")("triggerTimeOut_2p").getIntValue();
	ctrl()->writeRegister("triggerWindowEnd_2p", timeOut2p_10ns);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets modes of CoBo primary and secondary trigger logics.
 */
void CoBoNode::setTriggerModes(CCfg::CConfig & nodeCfg)
{
	setTriggerMode(nodeCfg);
	setSecondaryTriggerMode(nodeCfg);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets CoBo trigger mode.
 */
void CoBoNode::setTriggerMode(CCfg::CConfig & nodeCfg)
{
	// First, all trigger modes are disabled.
	ctrl()->writeField("triggerMask", "external", 0);
	ctrl()->writeField("triggerMask", "multiplicity", 0);
	ctrl()->writeField("triggerMask", "timing", 0);
	ctrl()->writeField("triggerMask", "manual", 0);
	ctrl()->writeField("triggerMask", "slowControl", 0);
	ctrl()->writeField("triggerMask", "external_2p", 0);
	ctrl()->writeField("triggerMask", "multiplicity_2p", 0);
	ctrl()->writeField("triggerMask", "timing_2p", 0);
	ctrl()->writeField("triggerMask", "manual_2p", 0);
	ctrl()->writeField("triggerMask", "slowControl_2p", 0);

	// The requested trigger mode is enabled.
	const std::string mode = nodeCfg("Module")("triggerMode").getStringValue();
	NODE_INFO() << "Selecting '" << mode << "' primary trigger mode";
	if (mode == "onRequest")
	{
		ctrl()->writeField("triggerMask", "manual", 1);
	}
	else if (mode == "onMultiplicity")
	{
		ctrl()->writeField("triggerMask", "multiplicity", 1);
	}
	else if (mode == "onSlowControl")
	{
		ctrl()->writeField("triggerMask", "slowControl", 1);
	}
	else if (mode == "periodically")
	{
		ctrl()->writeField("triggerMask", "timing", 1);
	}
	else if (mode == "external")
	{
		ctrl()->writeField("triggerMask", "external", 1);
	}
	else if (mode == "onRequest_2p")
	{
		ctrl()->writeField("triggerMask", "manual", 1);
		ctrl()->writeField("triggerMask", "manual_2p", 1);
	}
	else if (mode == "onMultiplicity_2p")
	{
		ctrl()->writeField("triggerMask", "multiplicity", 1);
		ctrl()->writeField("triggerMask", "multiplicity_2p", 1);
	}
	else if (mode == "onSlowControl_2p")
	{

		ctrl()->writeField("triggerMask", "slowControl_2p", 1);
		ctrl()->writeField("triggerMask", "slowControl", 1);
	}
	else if (mode == "periodically_2p")
	{
		ctrl()->writeField("triggerMask", "timing", 1);
		ctrl()->writeField("triggerMask", "timing_2p", 1);
	}
	else if (mode == "external_2p")
	{
		ctrl()->writeField("triggerMask", "external", 1);
		ctrl()->writeField("triggerMask", "external_2p", 1);
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 *  Sets mode of the CoBo secondary trigger.
 *  This secondary trigger logic block, when it triggers, sends a pulse on connector J6, pin 4 but does not trigger the event readout.
 *  The trigger block uses the same parameters as the primary trigger, but has a different trigger mask.
 */
void CoBoNode::setSecondaryTriggerMode(CCfg::CConfig & nodeCfg)
{
	// First, all trigger modes are disabled.
	try
	{
		ctrl()->writeField("triggerMask", "external_out", 0);
		ctrl()->writeField("triggerMask", "multiplicity_out", 0);
		ctrl()->writeField("triggerMask", "timing_out", 0);
		ctrl()->writeField("triggerMask", "manual_out", 0);
		ctrl()->writeField("triggerMask", "slowControl_out", 0);
	}
	catch (const CmdException & e)
	{
		NODE_WARN() << "Error setting secondary trigger mode: " << e.reason;
		return;
	}

	// The requested trigger mode is enabled.
	std::string mode = "noTrigger";
	try
	{
		try
		{
			mode = nodeCfg("Module")("secondaryTriggerMode").getStringValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			mode = nodeCfg("Module")("externalTriggerMode").getStringValue(); // old short-lived parameter name
		}
	}
	catch (const CCfg::CConfig::NotFound&)
	{
		mode = nodeCfg("Module")("triggerMode").getStringValue();
		NODE_WARN() << "Using same mode for primary and secondary trigger logics.";
	}

	NODE_DEBUG() << "Selecting '" << mode << "' secondary trigger mode";

	if (mode == "onRequest")
	{
		ctrl()->writeField("triggerMask", "manual_out", 1);
	}
	else if (mode == "onMultiplicity")
	{
		ctrl()->writeField("triggerMask", "multiplicity_out", 1);
	}
	else if (mode == "onSlowControl")
	{
		ctrl()->writeField("triggerMask", "slowControl_out", 1);
	}
	else if (mode == "periodically")
	{
		ctrl()->writeField("triggerMask", "timing_out", 1);
	}
	else if (mode == "external")
	{
		ctrl()->writeField("triggerMask", "external_out", 1);
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the delay from rising edge of SCW signal to start of multiplicity sliding window.
 * @param nodeCfg Configuration to extract the delay value from.
 */
void CoBoNode::setScwMultDelay(CCfg::CConfig & nodeCfg)
{
	uint64_t delay = 200;
	try
	{
		// Get delay value in units of CKR clock cycles.
		delay = nodeCfg("Module")("scwMultDelay").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
	}
	NODE_DEBUG() << "Setting 'ScwMultDelay' to " << delay << " CKR clock cycle(s)";
	ctrl()->writeRegister("ScwMultDelay", delay);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the delay between SCW low and SCR high.
 * @param nodeCfg Configuration to extract the delay value from.
 */
void CoBoNode::setScwScrDelay(CCfg::CConfig & nodeCfg)
{
	uint64_t delay_10ns = 0;
	try
	{
		// Get delay value in units of 10 nanoseconds.
		delay_10ns = nodeCfg("Module")("scwScrDelay").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound&)
	{
		try
		{
			delay_10ns = nodeCfg("Module")("writeReadDelay").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound&)
		{
			// Keep default value
		}
	}
	NODE_DEBUG() << "Setting 'ScwScrDelay' to " << delay_10ns*10 << " ns";
	ctrl()->writeRegister("ScwScrDelay", delay_10ns);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configures the CoBo circular buffers.
 */
void CoBoNode::setCircularBuffers(CCfg::CConfig & nodeCfg)
{
	uint32_t bufferStart[4], bufferEnd[4];
	// Get total size
	const uint32_t startOfMemory = getCircularBufferParameter(nodeCfg, "startOfMemory");
	const uint32_t endOfMemory = getCircularBufferParameter(nodeCfg, "endOfMemory");

	// Check if the location and size of each individual buffer is provided in configuration
	try
	{
		for (size_t asadIdx=0; asadIdx <4u; ++asadIdx)
		{
			bufferStart[asadIdx] = getAsAdCircularBufferParameter(nodeCfg, asadIdx, "startOfMemory");
			bufferEnd[asadIdx] = getAsAdCircularBufferParameter(nodeCfg, asadIdx, "endOfMemory");
		}
	}
	// Otherwise, compute size of each AsAd board's circular buffer depending on the number of boards
	catch (const CCfg::CConfig::NotFound &)
	{
		// Get number of active boards
		size_t numAsAdBoards = 0;
		for (size_t asadIdx=0; asadIdx <4u; ++asadIdx)
		{
			if (isAsadActive(asadIdx)) numAsAdBoards++;
		}
		if (3 == numAsAdBoards) numAsAdBoards  = 4;
		// Compute buffer sizes
		uint32_t totalSize = endOfMemory - startOfMemory;
		uint32_t bufferSize = totalSize/numAsAdBoards;
		bufferStart[0] = startOfMemory;
		for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
		{
			uint32_t currentSize = (isAsadActive(asadIdx))?bufferSize:0;
			bufferEnd[asadIdx] = bufferStart[asadIdx] + currentSize;
			if (asadIdx < 3) bufferStart[asadIdx+1] = bufferEnd[asadIdx];
		}
	}

	// Set circular buffer addresses
	try
	{
		for (size_t asadIdx = 0; asadIdx < 4u; ++asadIdx)
		{
			std::string regName = std::string("memStart_asad") + boost::lexical_cast<std::string>(asadIdx);
			ctrl()->writeRegister(regName, bufferStart[asadIdx]);
			regName = std::string("memEnd_asad") + boost::lexical_cast<std::string>(asadIdx);
			ctrl()->writeRegister(regName, bufferEnd[asadIdx]);
			regName = std::string("readPtr_asad") + boost::lexical_cast< std::string >(asadIdx);
			ctrl()->writeRegister(regName, bufferStart[asadIdx]);
		}
	}
	catch (const CmdException &)
	{
		ctrl()->writeRegister("memStart", startOfMemory);
		ctrl()->writeRegister("memEnd", endOfMemory);
		ctrl()->writeRegister("readPtr", startOfMemory);
	}

	// Set interrupt rate
	uint32_t eventsPerInterrupt = 1u;
	try
	{
		eventsPerInterrupt = getCircularBufferParameter(nodeCfg, "interruptRate");
	}
	catch (const CCfg::CConfig::NotFound &) {}
	try
	{
		ctrl()->writeRegister("interruptRate", eventsPerInterrupt);
	}
	catch (const CmdException &) {}

	// Set levels
	ctrl()->writeRegister("aFullLevel", getCircularBufferParameter(nodeCfg, "levelAlmostFullMemory"));
	ctrl()->writeRegister("aEmptyLevel", getCircularBufferParameter(nodeCfg, "levelAlmostEmptyMemory"));
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Selects which of the circular buffers should be read.
 * @param mask Mask of enabled buffers (1 bit per buffer).
 */
void CoBoNode::setCircularBuffersEnabled(uint16_t mask)
{
	backend->selectHwNodeByEndpoint(hwNode().endpoint);
	backend->setCircularBuffersEnabled(mask);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Configures the AsAd TCM device.
 *   - selects step voltage or current pulse mode: T/C (aka Calib)
 *   - sets capacitance: RNG (aka Rg): 1pF (0) or 11 pF (1) capacitor.
 *   - selects trigger mode (external or slow control): TRG
 *   - sets triggering delay: DLY
 *   - sets default output voltage
 *   - sets output voltage for external trigger mode
 *
 * @see Section 1.4 of GET-AS-002-0005: "AsAd Test and Calibration Management".
 * @param nodeCfg CoBo board configuration
 * @param asadId AsAd board identifier
 */
void CoBoNode::configureAsAdTcmDevice(CCfg::CConfig & nodeCfg, size_t asadId)
{
	// Check mode: either the pulser signal goes directly to the AGET test input (T) or it goes through a calibrated external capacitor (C).
	const bool TC = asadPulserCalibMode(nodeCfg, asadId);
	// Check which of the 1 pF and 11 pF external capacitances to use.
	const bool RNG = asadPulser11pFCapacitorSelected(nodeCfg, asadId);
	// Check whether pulser should be triggered by slow control modification of TCM2 or by external signal from ISP1
	const bool TRG = asadPulserTriggerMode(nodeCfg, asadId);
	// Check triggering delay.
	const uint8_t DLY = asadPulserTriggerDelay(nodeCfg, asadId);

	NODE_DEBUG() << "Configuring generator of AsAd no. " << asadId << " (T/C=" << TC << ", RNG=" << RNG << " TRG=" << TRG << " DLY=" << (uint16_t) DLY << ')';
	asad()->writeField("TCM0", "Calib", TC);
	asad()->writeField("TCM0", "Rg", RNG);
	asad()->writeField("TCM0", "Delay", DLY);

	// Enable AsAd TCM device slow control triggering mode
	asad()->writeField("TCM0", "Trg", 0);

	static const float Vref_mV = 1102.0f;
	// Check default voltage
	const float defaultVoltage_mV = asadPulserDefaultVoltage_mV(nodeCfg);
	// Check pulse voltage
	const float voltage_mV = asadPulserVoltage_mV(nodeCfg);

	// Set default output voltage
	NODE_DEBUG() << "Setting AsAd generator default output voltage to " <<  (Vref_mV + defaultVoltage_mV) << " mV";
	asadPulserMgr()->setDefaultPulserVoltage(defaultVoltage_mV);

	// Enable AsAd TCM device external triggering mode if needed
	if (TRG)
	{
		NODE_DEBUG() << "Setting AsAd generator output voltage to " <<  (Vref_mV + voltage_mV) << " mV";
		asad()->writeField("TCM0", "Trg", TRG);
		asadPulserMgr()->triggerPulser(voltage_mV);
	}
	try
	{
		asadPulserMgr()->configureExternalPulser(TRG, voltage_mV);
	}
	catch (const ::Ice::OperationNotExistException & e) // This operation was added later on
	{
		NODE_WARN() << e.what();
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Stores the configuration parameters for the periodic AsAd pulser, for later use.
 * Configures periodic AsAd pulser.
 */
void CoBoNode::setAsAdPeriodicPulser(CCfg::CConfig & nodeCfg)
{
	// Store parameters
	try
	{
		// Enabled?
		configData_.asadPulserPeriodicModeEnabled = nodeCfg("Generator")("periodicModeEnabled").getBoolValue();
		// Period
		configData_.asadPulserPeriod_ms = nodeCfg("Generator")("period_ms").getIntValue();
		// Amplitude
		configData_.asadPulserVoltage_mV = asadPulserVoltage_mV(nodeCfg);
		// 2p mode
		try
		{
			configData_.asadPulserDoubleModeEnabled = nodeCfg("Generator")("doublePulseEnabled").getBoolValue();
			configData_.asadPulserVoltage2_mV = nodeCfg("Generator")("secondPulseAmplitude").getRealValue();
		}
		catch (CCfg::CConfig::NotFound & e)
		{
			NODE_WARN() << e.what();
			configData_.asadPulserDoubleModeEnabled = false;
		}
		// Pseudo-random mode
		try
		{
			configData_.asadPulserRandomModeEnabled = nodeCfg("Generator")("poissonModeEnabled").getBoolValue();
		}
		catch (CCfg::CConfig::NotFound & e)
		{
			NODE_WARN() << e.what();
			configData_.asadPulserRandomModeEnabled = false;
		}
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		NODE_WARN() << e.what();
		configData_.asadPulserPeriodicModeEnabled = false;
		configData_.asadPulserPeriod_ms = 1000u;
		configData_.asadPulserVoltage_mV = -1000;
		configData_.asadPulserDoubleModeEnabled = false;
	}

	// Configure pulser manager
	if (configData_.asadPulserPeriodicModeEnabled)
	{
		asadPulserMgr()->configurePeriodicPulser(configData_.asadPulserVoltage_mV, configData_.asadPulserPeriod_ms);
		try
		{
			asadPulserMgr()->configureDoublePulserMode(configData_.asadPulserDoubleModeEnabled?configData_.asadPulserVoltage2_mV:0);
		}
		catch (const ::Ice::OperationNotExistException & e) // This operation was added later on
		{
			if (configData_.asadPulserDoubleModeEnabled)
			{
				throw mdaq::utl::CmdException(e.what());
			}
		}
		try
		{
			asadPulserMgr()->setRandomPulserEnabled(configData_.asadPulserRandomModeEnabled);
		}
		catch (const ::Ice::OperationNotExistException & e) // This operation was added later on
		{
			if (configData_.asadPulserRandomModeEnabled)
			{
				throw mdaq::utl::CmdException(e.what());
			}
		}
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the threshold and the size of the sliding window for the global multiplicity.
 */
void CoBoNode::setMultiplicityMode(CCfg::CConfig nodeCfg)
{
	// Set the multiplicity threshold
	const uint32_t threshold = nodeCfg("Module")("multiplicityThreshold").getIntValue();
	NODE_DEBUG() << "Setting multiplicity threshold to " << threshold;
	ctrl()->writeRegister("multThreshold", threshold);
	// Set the number of samples to sum in the sliding window
	const uint32_t windowWidth = nodeCfg("Module")("multWindowSize").getIntValue();
	NODE_DEBUG() << "Setting width of multiplicity sliding window to " << windowWidth;
	ctrl()->writeRegister("multWindSize", windowWidth);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the multiplicity threshold for the second trigger.
 */
void CoBoNode::setMultiplicityThreshold2p(CCfg::CConfig nodeCfg)
{
	const uint32_t secondThreshold = nodeCfg("Module")("multiplicityThreshold_2p").getIntValue();
	NODE_DEBUG() << "Setting 2p mode second multiplicity threshold to " << secondThreshold;
	ctrl()->writeRegister("multThreshold_2p", secondThreshold);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 *   Sets AsAd monitoring VDD limits.
 */
void CoBoNode::setAsAdVddLimits(CCfg::CConfig nodeCfg, size_t asadId)
{
	// Get limits in V
	double lowerVoltage = 0, upperVoltage = 7.0;
	try
	{
		lowerVoltage = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsVDD")("lowerValue").getRealValue();
		upperVoltage = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsVDD")("upperValue").getRealValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		lowerVoltage = asadConfig(nodeCfg, asadId)("Monitoring")("VddLowLimit").getRealValue();
		upperVoltage = asadConfig(nodeCfg, asadId)("Monitoring")("VddHighLimit").getRealValue();
	}
	NODE_DEBUG() << "Setting monitoring range of VDD to ["
			<< lowerVoltage << ", " << upperVoltage << "] V for AsAd no. " << asadId;

	// Convert to straight binary
	const uint8_t lowerLimit = 512. * lowerVoltage / 7.524;
	const uint8_t upperLimit = 512. * upperVoltage / 7.524;

	// Set registers
	asad()->writeRegister("VDDLowLimit", lowerLimit);
	asad()->writeRegister("VDDHighLimit", upperLimit);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 *   Sets AsAd monitoring VDD1 limits.
 */
void CoBoNode::setAsAdVdd1Limits(CCfg::CConfig nodeCfg, size_t asadId)
{
	// Get limits in V
	double lowerVoltage = 0, upperVoltage = 7.0;
	try
	{
		lowerVoltage = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsVDD1")("lowerValue").getRealValue();
		upperVoltage = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsVDD1")("upperValue").getRealValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		lowerVoltage = asadConfig(nodeCfg, asadId)("Monitoring")("Vdd1LowLimit").getRealValue();
		upperVoltage = asadConfig(nodeCfg, asadId)("Monitoring")("Vdd1HighLimit").getRealValue();
	}
	NODE_DEBUG() << "Setting monitoring range of VDD1 to ["
			<< lowerVoltage << ", " << upperVoltage << "] V for AsAd no. " << asadId;

	// Convert to straight binary
	const uint8_t lowerLimit = 256. * (lowerVoltage - 0.06) / 7.;
	const uint8_t upperLimit = 256. * (upperVoltage - 0.06) / 7.;

	// Set registers
	asad()->writeRegister("VDD1LowLimit", lowerLimit);
	asad()->writeRegister("VDD1HighLimit", upperLimit);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 *   Sets AsAd monitoring IDD limits.
 */
void CoBoNode::setAsAdIddLimits(CCfg::CConfig nodeCfg, size_t asadId)
{
	// Get limits in A
	double lowerCurrent = 0, upperCurrent = 3.0;
	try
	{
		lowerCurrent = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsIDD")("lowerValue").getRealValue();
		upperCurrent = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsIDD")("upperValue").getRealValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		lowerCurrent = asadConfig(nodeCfg, asadId)("Monitoring")("IddLowLimit").getRealValue();
		upperCurrent = asadConfig(nodeCfg, asadId)("Monitoring")("IddHighLimit").getRealValue();
	}
	NODE_DEBUG() << "Setting monitoring range of IDD to ["
			<< lowerCurrent << ", " << upperCurrent << "] A for AsAd no. " << asadId;

	// Convert to straight binary
	const uint8_t lowerLimit = 193.28 * lowerCurrent / 2.28;
	const uint8_t upperLimit = 193.28 * upperCurrent / 2.28;

	// Set registers
	asad()->writeRegister("IDDLowLimit", lowerLimit);
	asad()->writeRegister("IDDHighLimit", upperLimit);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 *   Set AsAd monitoring internal temperature limits.
 */
void CoBoNode::setAsAdIntTempLimits(CCfg::CConfig nodeCfg, size_t asadId)
{
	// Get limits in °C [-128:127]
	uint8_t lowerLimit = 0, upperLimit = 50;
	try
	{
		lowerLimit = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsTint")("lowerValue").getIntValue();
		upperLimit = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsTint")("upperValue").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		lowerLimit = asadConfig(nodeCfg, asadId)("Monitoring")("InternalTempLowLimit").getIntValue();
		upperLimit = asadConfig(nodeCfg, asadId)("Monitoring")("InternalTempHighLimit").getIntValue();
	}
	NODE_DEBUG() << "Setting monitoring range of temperature far from ADC to ["
			<< (uint16_t) lowerLimit << ", " << (uint16_t) upperLimit << "] deg C for AsAd no. " << asadId;
	asad()->writeRegister("intTempLowLimit", lowerLimit);
	asad()->writeRegister("intTempHighLimit", upperLimit);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 *   Set AsAd monitoring external temperature limits.
 */
void CoBoNode::setAsAdExtTempLimits(CCfg::CConfig nodeCfg, size_t asadId)
{
	// Get limits in °C [-128:127]
	uint8_t lowerLimit = 0, upperLimit = 50;
	try
	{
		lowerLimit = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsText")("lowerValue").getIntValue();
		upperLimit = asadConfig(nodeCfg, asadId)("Monitoring")("LimitsText")("upperValue").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		lowerLimit = asadConfig(nodeCfg, asadId)("Monitoring")("ExternalTempLowLimit").getIntValue();
		upperLimit = asadConfig(nodeCfg, asadId)("Monitoring")("ExternalTempHighLimit").getIntValue();
	}
	NODE_DEBUG() << "Setting monitoring range of temperature close to ADC to ["
			<< (uint16_t) lowerLimit << ", " << (uint16_t) upperLimit << "] deg C for AsAd no. " << asadId;
	asad()->writeRegister("extTempLowLimit", lowerLimit);
	asad()->writeRegister("extTempHighLimit", upperLimit);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Reads AsAd VDD voltage (see AsAd Power Supply and Monitoring document).
 * @param asadId Index of the AsAd board (between 0 and 3).
 * @return Value of VDD [V].
 */
double CoBoNode::readAsAdVdd(size_t asadId)
{
	uint64_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		// Select AsAd board
		ctrl()->writeField("asadSelect", "d", asadId);
		// Read register
		binValue = (uint64_t) asad()->readRegister("VDD");
	}
	// Convert from straight binary to volts
	double value = (binValue != 0) ? (2.28 / 256) * (3.3 / 2) * binValue : 0.0;
	NODE_DEBUG() << "Read AsAd " << asadId << " VDD voltage: " << value << " V";
	return value;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Reads the VDD1 register of an AsAd board (see AsAd Power Supply and Monitoring document).
 * @param asadId Index of the AsAd board (between 0 and 3).
 * @return Value of VDD1 [V].
 */
double CoBoNode::readAsAdVdd1(size_t asadId)
{
	uint64_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		// Select AsAd board
		ctrl()->writeField("asadSelect", "d", asadId);
		// Read register
		binValue = (uint64_t) asad()->readRegister("VDD1");
	}
	// Convert to volts
	double value = (binValue != 0) ? (60e-3 + 7 * ((double) binValue) / 256) : 0.0;
	NODE_DEBUG() << "Read AsAd " << asadId << " VDD1 voltage: " << value << " V";
	return value;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Reads the IDD register of an AsAd board (see AsAd Power Supply and Monitoring document).
 * @param asadId Index of the AsAd board (between 0 and 3).
 * @return Value of IDD [A].
 */
double CoBoNode::readAsAdIdd(size_t asadId)
{
	uint64_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		// Select AsAd board
		ctrl()->writeField("asadSelect", "d", asadId);
		// Read register
		binValue = (uint64_t) asad()->readRegister("IDD");
	}
	// Convert to A
	double value = binValue != 0 ? (2.28 / 256) * (1 / 0.755) * binValue : 0.0;
	NODE_INFO() << "Read AsAd " << asadId << " IDD current: " << value << " A";
	return value;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Reads AsAd on-chip temperature (see AsAd Power Supply and Monitoring document).
 * @param asadId Index of the AsAd board (between 0 and 3).
 * @return Temperature [° C].
 */
int8_t CoBoNode::readAsAdIntTemp(size_t asadId)
{
	uint64_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		// Select AsAd board
		ctrl()->writeField("asadSelect", "d", asadId);
		// Read register
		binValue = (uint64_t) asad()->readRegister("internalTemperature");
	}
	// Convert to degrees C
	uint8_t value = ((binValue & 0x80) == 0x8)? ~(binValue - 1) : binValue;
	return value;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Reads  AsAd ADC-chip surroundings temperature (see AsAd Power Supply and Monitoring document).
 * @param asadId Index of the AsAd board (between 0 and 3).
 * @return Temperature [° C].
 */
int8_t CoBoNode::readAsAdExtTemp(size_t asadId)
{
	uint64_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		// Select AsAd board
		ctrl()->writeField("asadSelect", "d", asadId);
		// Read register
		binValue = (uint64_t) asad()->readRegister("externalTemperature");
	}
	// Convert to degrees C
	uint8_t value = ((binValue & 0x80) == 0x8)? ~(binValue - 1) : binValue;
	return value;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the ISPM multiplexer to debug the signals on AsAd.
 * See Fig. 14 and Tables 14-20 page 17 of "AsAd Firmware Description".
 * See "AsAd External Instruments Interface".
 * @param nodeCfg Configuration to retrieve the parameters from.
 * @param asadId Index of the AsAd board (between 0 and 3).
 */
void CoBoNode::setAsAdInspectionManager(CCfg::CConfig nodeCfg, const std::size_t & asadId)
{
	// Lookup register value.
	std::string device, isp1, isp2;
	CoBoNode::getAsadInspectionLinesParameters(nodeCfg, asadId, device, isp1, isp2);
	const uint64_t value = get::asad::utl::buildIspmRegisterValue(device, isp1, isp2);

	// Set register
	NODE_DEBUG() << "Setting AsAd no. " << asadId << " inspection manager to " << std::showbase << std::hex << value << std::dec <<
			" for signals (" << isp1 << ", " << isp2 << ") of device " << device;
	asad()->writeRegister("ISPM", value);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Checks desired readout mode in configuration.
 * Returns true (full readout mode) if all active chips have selected this mode.
 */
bool CoBoNode::isFullReadoutModeEnabled(CCfg::CConfig & nodeCfg)
{
	// Check readout mode of all chips
	size_t partialReadoutCount = 0;
	size_t activeChipCount = 0;
	for (size_t asadIdx=0; asadIdx < 4; ++asadIdx)
	{
		for (size_t agetIdx=0; agetIdx < 4; ++agetIdx)
		{
			if (not isAgetActive(asadIdx, agetIdx)) continue;

			activeChipCount++;
			if (not agetConfig(nodeCfg, asadIdx, agetIdx)("Global")("Reg1")("isAllChannelRead").getBoolValue())
				partialReadoutCount++;
		}
	}

	// Ensure all chips are in the same readout mode
	bool fullReadoutModeEnabled = true;
	if (partialReadoutCount > 0)
	{
		if (partialReadoutCount == activeChipCount)
		{
			fullReadoutModeEnabled = false;
		}
		else
		{
			std::ostringstream msg;
			msg << "All AGET chips of CoBo board '" << id().toString()
				<< "' should be in the same (full or partial) readout mode! Check the 'isAllChannelRead' parameters.";
			throw mdaq::utl::CmdException(msg.str());
		}
	}
	return fullReadoutModeEnabled;
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets (partial/full) readout mode CoBo register.
 * All AGET chips within a CoBo should be configured accordingly (see AGET register 1 bit 14).
 */
void CoBoNode::setReadoutMode(CCfg::CConfig & nodeCfg)
{
	bool fullReadoutMode = isFullReadoutModeEnabled(nodeCfg);

	NODE_INFO() << "Selecting " << (fullReadoutMode?"full":"partial") << " readout mode";
	ctrl()->writeField("pipeCtrl", "partialReadOut", fullReadoutMode?0:1);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the readout depth (number of AGET time buckets to read).
 *
 * If the corresponding parameter is not found, the number of buckets to read is set to 512 or 256 instead.
 * FIXME: check for overflow taking readout offset into account
 */
void CoBoNode::setReadoutDepth(CCfg::CConfig & nodeCfg)
{
	const bool is2pModeEnabled = nodeCfg("Module")("enableMem2pMode").getBoolValue();
	uint16_t readoutDepth = is2pModeEnabled ? 256 : 512 ;
	try
	{
		readoutDepth = nodeCfg("Module")("readoutDepth").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		NODE_WARN() << e.what() << " Using a readout depth of " << readoutDepth << " bucket(s) instead.";
	}
	if (is2pModeEnabled and readoutDepth > 256)
	{
		NODE_ERROR() << readoutDepth << " is not a valid readout depth in 2p mode. Using 256 bucket(s) instead.";
		readoutDepth = 256;
	}
	else if (not is2pModeEnabled and readoutDepth > 512)
	{
		NODE_ERROR() << readoutDepth << " is not a valid readout depth. Using 512 bucket(s) instead.";
		readoutDepth = 512;
	}
	NODE_DEBUG() << "Setting readout depth to " << readoutDepth << " time bucket(s)";
	ctrl()->writeRegister("readDepth", readoutDepth);
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Sets the 'readOffset' CoBo register indicating to CoBo the offset used by the AGET chips
 * for the readout pointer (AGET register 12) and to be copied in the 'readOffset' field
 * of the data frame header.
 *
 * Sets the CoBo 'readDepth' register to the value given to the register 12 of the first active AGET chip.
 */
void CoBoNode::setReadoutOffset(CCfg::CConfig & nodeCfg)
{
	// Find first active chip
	const size_t fullAgetIdx = firstActiveChip();
	const size_t asadIdx = fullAgetIdx/4;
	const size_t agetIdx = fullAgetIdx%4;

	// Get offset value
	const uint16_t readoutOffset = agetConfig(nodeCfg, asadIdx, agetIdx)("Global")("ReadoutPointerOffset").getIntValue();
	MSG_IF(readoutOffset>0u) << "Setting readout offset to " << readoutOffset << " time bucket(s)";

	// Set register
	try
	{
		ctrl()->writeField("readOffset", "value", readoutOffset);
	}
	catch (const CmdException & e)
	{
		NODE_WARN() << e.reason; // This register was only described on Nov 25th, 2013
	}
}
// --------------------------------------------------------------------------------------------------------------------
/**
 * Initializes 2p mode readout.
 */
void CoBoNode::aget2pInitNow()
{
	NODE_DEBUG() << "Initializing 2p mode readout";
	ctrl()->writeRegister("aget_2p_initNow", 1);
	sleep_ms(20);
	ctrl()->writeRegister("aget_2p_initNow", 0);
}
// --------------------------------------------------------------------------------------------------------------------
CoBoNode::ConfigData::ConfigData()
	: asadPulserPeriodicModeEnabled(false),
	  asadPulserPeriod_ms(1000u),
	  asadPulserVoltage_mV(-1000),
	  asadPulserDoubleModeEnabled(false),
	  asadPulserRandomModeEnabled(false),
	  asadPulserVoltage2_mV(-600)
{
};
// --------------------------------------------------------------------------------------------------------------------
} /* namespace rc */
} /* namespace get */
