/* =====================================================================================================================
 * GetEccImpl.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class get::GetEccImpl
 * Created on: Jan 11, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
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
#include "get/GetEccImpl.h"
#include "get/rc/BackEnd.h"
#include "get/rc/SubConfigIdSet.h"
#include "get/rc/ConfigIdSet.h"
#include "mdaq/cfg/NodeAtTop.h"
#include "mdaq/DefaultPortNums.h"
#include "utl/Logging.h"
#include ICE_HEADER(ObjectAdapter.h)
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace get
{
//_________________________________________________________________________________________________
GetEccImpl::GetEccImpl() : ::mdaq::EccImpl()
{
}
//_________________________________________________________________________________________________
GetEccImpl::~GetEccImpl()
{
}
//_________________________________________________________________________________________________
rc::BackEnd & GetEccImpl::backEnd()
{
	if (backEndPtr.get() == 0)
	{
		backEndPtr.reset(new rc::BackEnd());
	}
	return dynamic_cast< rc::BackEnd & >(*backEndPtr.get());
}
//_________________________________________________________________________________________________
/** Sets up the manger of configurations available to the state machine.
 * @param repoUrl Path of the directory or URL of the database containing the configurations.
 */
void GetEccImpl::setConfigManager(const std::string & repoUrl)
{
	backEnd().setConfigManager(repoUrl);
}
//_________________________________________________________________________________________________
void GetEccImpl::createAlarmLogger(const ::utl::net::SocketAddress & address)
{
	backEnd().createAlarmLogger(address);
}
//_________________________________________________________________________________________________
void GetEccImpl::daqStart(const Ice::Current &)
{
	backEnd().daqStart();
}
//_________________________________________________________________________________________________
void GetEccImpl::daqStop( const Ice::Current &)
{
	backEnd().daqStop();
}
//_________________________________________________________________________________________________
void GetEccImpl::describeNode(const std::string & nodeType, const std::string & instanceId, const std::string & configId, const Ice::Current &)
{
	backEnd().describeNode(nodeType, instanceId, configId);
}
//_________________________________________________________________________________________________
void GetEccImpl::prepareNode(const std::string & nodeType, const std::string & instanceId, const std::string & configId, const Ice::Current &)
{
	backEnd().prepareNode(nodeType, instanceId, configId);
}
//_________________________________________________________________________________________________
void GetEccImpl::configureNode(const std::string & nodeType, const std::string & instanceId, const std::string & configId, const Ice::Current &)
{
	backEnd().configureNode(nodeType, instanceId, configId);
}
//_________________________________________________________________________________________________
/*=============================================================
 * CoBo specific
 =============================================================*/
void GetEccImpl::systemReset(const Ice::Current& ) {}
void GetEccImpl::coboSetWriteEnable(const Ice::Current& ) {}
void GetEccImpl::coboSetReadEnable(const Ice::Current& ) {}
int  GetEccImpl::coboReadMultiplicity(const Ice::Current& ) {return 0;}
//_________________________________________________________________________________________________
/**
 * Sets CoBo writing clock CKW frequency.
 * @param frequency_MHz CKW requency [MHz].
 * @param pllConfigEnabled Debugging flag to disable configuration of the CoBo PLL device.
 * @param registers Map of non-default PLL register values.
 * @param c Information about the current method invocation.
 */
void GetEccImpl::coboSetWritingClockFrequency(float frequency_MHz, bool pllConfigEnabled, const cobo::PllRegisterMap & registers, const Ice::Current & /* c */)
{
	backEnd().coboSetWritingClockFrequency(frequency_MHz, pllConfigEnabled, registers);
}
//_________________________________________________________________________________________________
/**
 * Sets CoBo reading clock CKR frequency.
 * @param frequency_MHz CKR frequency [MHz].
 * @param c Information about the current method invocation.
 * The only frequency value allowed is 25 MHz.
 */
void GetEccImpl::coboSetReadingClockFrequency(float frequency_MHz, const Ice::Current &)
{
	backEnd().coboSetReadingClockFrequency(frequency_MHz);
}
//_________________________________________________________________________________________________
/**
 * Switches current CoBo node to slow or fast control mode.
 * @param slowMode Whether to switch from fast to slow (true) or from slow to fast control mode.
 * @param c Information about the current method invocation.
 */
void GetEccImpl::coboSwitchToSlowControlMode(bool slowMode, const Ice::Current & /* c */)
{
	backEnd().coboSwitchToSlowControlMode(slowMode);
}
//_________________________________________________________________________________________________
/**
 * Selects single channel (per chip) of current CoBo node to always read.
 * @param chanId Index of the single channel to read (0 to 67) or 68 to read none.
 * @param c Information about the current method invocation.
 */
void GetEccImpl::coboSelectSingleAgetReadoutChannel(::Ice::Long chanId, const Ice::Current & /* c */)
{
	backEnd().coboSelectSingleReadoutChannel(chanId);
}
//_________________________________________________________________________________________________
/**
 * Selects single channel (per chip) of current CoBo node for test mode.
 * @param chanId Index of the single test channel to read (0 to 67). If 68, select none.
 * @param c Information about the current method invocation.
 */
void GetEccImpl::coboSelectSingleAgetTestChannel(::Ice::Long chanId, const Ice::Current & /* c */)
{
	backEnd().coboSelectSingleTestChannel(chanId);
}
//_________________________________________________________________________________________________
/**
 * Modifies LSB threshold of a single channel (per chip) of current CoBo node for test mode.
 * @param chanId Index of the single test channel to read (0 to 67). If 68, set all.
 * @param threshold Least significant bits of channel threshold (from 0 to 15).
 * @param c Information about the current method invocation.
 */
void GetEccImpl::coboSetAgetChannelLsbThreshold(::Ice::Long chanId, ::Ice::Long threshold, const Ice::Current & /* c */)
{
	backEnd().coboSetAgetChannelLsbThreshold(chanId, threshold);
}
//_________________________________________________________________________________________________
/**
 * @deprecated Use triggerAsAdPulser instead.
 */
void GetEccImpl::asadInjectCharge(Ice::Long voltage_mV, const Ice::Current & c)
{
	triggerAsAdPulser(voltage_mV, c);
}
//_________________________________________________________________________________________________
void GetEccImpl::setAsAdPulserDefaultVoltage(Ice::Long voltage_mV, const Ice::Current &)
{
	backEnd().asadPulserMgr()->setDefaultPulserVoltage(voltage_mV);
}
//_________________________________________________________________________________________________
void GetEccImpl::triggerAsAdPulser(Ice::Long voltage_mV, const Ice::Current &)
{
	backEnd().asadPulserMgr()->triggerPulser(voltage_mV);
}
//_________________________________________________________________________________________________
void GetEccImpl::configureAsAdExternalPulser(bool enabled, Ice::Long voltage_mV, const Ice::Current &)
{
	backEnd().asadPulserMgr()->configureExternalPulser(enabled, voltage_mV);
}
//_________________________________________________________________________________________________
void GetEccImpl::configureAsAdPeriodicPulser(Ice::Long voltage_mV, Ice::Long period_ms, const Ice::Current &)
{
	backEnd().asadPulserMgr()->configurePeriodicPulser(voltage_mV, period_ms);
}
//_________________________________________________________________________________________________
void GetEccImpl::setRandomAsAdPulserEnabled(bool poissonMode, const Ice::Current &)
{
	try
	{
		backEnd().asadPulserMgr()->setRandomPulserEnabled(poissonMode);
	}
	catch (const Ice::LocalException & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
}
//_________________________________________________________________________________________________
void GetEccImpl::configureAsAdDoublePulserMode(Ice::Long voltage_mV, const Ice::Current &)
{
	try
	{
		backEnd().asadPulserMgr()->configureDoublePulserMode(voltage_mV);
	}
	catch (const Ice::LocalException & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
}
//_________________________________________________________________________________________________
void GetEccImpl::startAsAdPeriodicPulser(const Ice::Current &)
{
	backEnd().asadPulserMgr()->startPeriodicPulser();
}
//_________________________________________________________________________________________________
void GetEccImpl::stopAsAdPeriodicPulser(const Ice::Current &)
{
	backEnd().asadPulserMgr()->stopPeriodicPulser();
}
//__________________________________________________________________________________________________
/**
 * Returns proxy to alarm service of current CoBo node.
 */
get::mt::AlarmServicePrx GetEccImpl::alarmService(const Ice::Current &)
{
	return backEnd().alarmService();
}
//__________________________________________________________________________________________________
/**
 * Subscribes ECC to alarms from current CoBo node.
 */
void GetEccImpl::subscribeAlarms(bool enabled, const Ice::Current &)
{
	backEnd().subscribeAlarms(enabled);
}
//_________________________________________________________________________________________________
void GetEccImpl::setAsAdAlarmMonitoringEnabled(bool enabled, const Ice::Current &)
{
	backEnd().asadAlarmMonitor()->setAsAdAlarmMonitoringEnabled(enabled);
}
//__________________________________________________________________________________________________
/**
 * Returns proxy to LED device manager of current CoBo node.
 */
cobo::LedManagerPrx GetEccImpl::ledManager(const Ice::Current &)
{
	return backEnd().ledManager();
}
//_____________________________________AGet specific_______________________________________________
void GetEccImpl::agetReadHitRegister(std::vector<long>&, const Ice::Current& ) {}
void GetEccImpl::agetWriteHitRegister(const ::mdaq::LongList&, const Ice::Current& ) {}
void GetEccImpl::agetCalibrateSlowControlDelay(const Ice::Current& ) {}
void GetEccImpl::agetCalibrateFastControlDelay(const Ice::Current& ) {}
//_________________________________________________________________________________________________
/**
 * Sets the ISPM multiplexer to debug the signals on AsAd.
 * See Fig. 14 and Tables 14-20 page 17 of "AsAd Firmware Description".
 * See "AsAd External Instruments Interface".
 * @param moduleName Name of device.
 * @param isp1SignalName Name of ISP1 signal.
 * @param isp2SignalName Name of ISP2 signal.
 */
void GetEccImpl::asadSetISPM(const std::string & moduleName, const std::string & isp1SignalName, const std::string & isp2SignalName, const Ice::Current &)
{
	backEnd().asadSetISPM(moduleName, isp1SignalName, isp2SignalName);
}
//__________________________________________________________________________________________________
/**
 * Sets up hardware topology.
 */
void GetEccImpl::describe(const std::string & configID, const Ice::Current &)
{
	try
	{
		backEnd().describe(configID);
	}
	catch (const get::rc::SM::Exception & e)
	{
		throw mdaq::utl::CmdException(e.errorMessage);
	}
}
//__________________________________________________________________________________________________
/**
 * Inter-board communications, clock synchronizations, calibrations
 */
void GetEccImpl::prepare(const std::string & configID, const Ice::Current &)
{
	try
	{
		backEnd().prepare(configID);
	}
	catch (const rc::SM::Exception & e)
	{
		throw mdaq::utl::CmdException(e.errorMessage);
	}
}
//__________________________________________________________________________________________________
/**
 * Resets state machine to its initial state.
 */
void GetEccImpl::reset(const Ice::Current &)
{
	try
	{
		backEnd().reset();
	}
	catch (const rc::SM::Exception & e)
	{
		throw mdaq::utl::CmdException(e.errorMessage);
	}
}
//__________________________________________________________________________________________________
/**
 * Loads physics parameters into the registers of the boards (thresholds, ...).
 */
void GetEccImpl::configure(const std::string & configID, const std::string & dataRouter, const Ice::Current &)
{
	try
	{
		backEnd().configure(configID, dataRouter);
	}
	catch (const rc::SM::Exception & e)
	{
		throw mdaq::utl::CmdException(e.errorMessage);
	}
}
//__________________________________________________________________________________________________
/**
 * Starts the run, i.e. allows hardware to deliver data.
 */
void GetEccImpl::start(const Ice::Current &)
{
	try
	{
		backEnd().start();
	}
	catch (const rc::SM::Exception & e)
	{
		throw mdaq::utl::CmdException(e.errorMessage);
	}
}
//__________________________________________________________________________________________________
/**
 * Ends the run, preventing hardware from delivering data.
 * Leads to Ready state.
 */
void GetEccImpl::stop(const Ice::Current &)
{
	try
	{
		backEnd().stop();
	}
	catch (const rc::SM::Exception & e)
	{
		throw mdaq::utl::CmdException(e.errorMessage);
	}
}
//__________________________________________________________________________________________________
/**
 * Pauses the run, i.e. forbids hardware to deliver data.
 * Leads to Paused state.
 */
void GetEccImpl::pause(const Ice::Current &)
{
	backEnd().pause();
}
//__________________________________________________________________________________________________
/**
 * Resumes the run, i.e. allows again hardware to deliver data.
 */
void GetEccImpl::resume(const Ice::Current &)
{
	backEnd().resume();
}
//__________________________________________________________________________________________________
/**
 * Invalidates physics parameters setup in the boards.
 * Leads to Idle state (so physics parameters will have to be setup again).
 */
void GetEccImpl::breakup(const Ice::Current &)
{
	backEnd().breakup();
}
//__________________________________________________________________________________________________
rc::Status toIceStatus(const ::get::rc::SM::Status & status)
{
	rc::Status value;
	switch (status.state)
	{
	case rc::SM::IDLE:
		value.s = get::rc::Idle;
		break;
	case rc::SM::DESCRIBED:
		value.s = get::rc::Described;
		break;
	case rc::SM::PREPARED:
		value.s = get::rc::Prepared;
		break;
	case rc::SM::READY:
		value.s = get::rc::Ready;
		break;
	case rc::SM::RUNNING:
		value.s = get::rc::Running;
		break;
	case rc::SM::PAUSED:
		value.s = get::rc::Paused;
		break;
	case rc::SM::OFF:
	default:
		value.s = get::rc::Off;
		break;
	}
	switch (status.error)
	{
	case rc::SM::WHEN_DESCRIBE:
		value.e = get::rc::WhenDescribe;
		break;
	case rc::SM::WHEN_PREPARE:
		value.e = get::rc::WhenPrepare;
		break;
	case rc::SM::WHEN_CONFIGURE:
		value.e = get::rc::WhenConfigure;
		break;
	case rc::SM::WHEN_START:
		value.e = get::rc::WhenStart;
		break;
	case rc::SM::WHEN_STOP:
		value.e = get::rc::WhenStop;
		break;
	case rc::SM::WHEN_PAUSE:
		value.e = get::rc::WhenPause;
		break;
	case rc::SM::WHEN_RESUME:
		value.e = get::rc::WhenResume;
		break;
	case rc::SM::WHEN_BREAKUP:
		value.e = get::rc::WhenBreakup;
		break;
	case rc::SM::WHEN_RESET:
		value.e = get::rc::WhenReset;
		break;
	case rc::SM::NO_ERR:
	default:
		value.e = get::rc::NoErr;
		break;
	}
	return value;
}
//__________________________________________________________________________________________________
void GetEccImpl::getStatus(rc::Status & s, const Ice::Current &)
{
	get::rc::SM::Status status;
	backEnd().getStateMachineStatus(status);
	s = toIceStatus(status);
}
//__________________________________________________________________________________________________
std::string GetEccImpl::getDataLinks(const Ice::Current &)
{
	return backEnd().dataLinks().toXml();
}
//_________________________________________________________________________________________________
std::string GetEccImpl::addSubConfig(const std::string & subConfigType, const std::string & subConfig, const std::string & subConfigId, const Ice::Current &)
{
	try
	{
		return backEnd().addSubConfig(subConfigType, subConfig, subConfigId);
	}
	catch (const std::exception & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
	return "";
}
//__________________________________________________________________________________________________
std::string GetEccImpl::getSubConfig(const std::string & subConfigType, const std::string & subConfigId, const Ice::Current &)
{
	try
	{
		return backEnd().getSubConfig(subConfigType, subConfigId);
	}
	catch (const std::exception & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
	return "";
}
//__________________________________________________________________________________________________
std::string GetEccImpl::getConfigSet(const std::string & xmlConfigId, const Ice::Current &)
{
	try
	{
		get::rc::ConfigId configId;
		configId.fromXml(xmlConfigId);
		return backEnd().getConfigSet(configId);
	}
	catch (const std::exception & e)
	{
		throw ::mdaq::utl::CmdException(e.what());
	}
	return "";
}
//__________________________________________________________________________________________________
/**
 * Returns all configuration IDs that are available to the ECC.
 * @return A sequence of string identifiers.
 */
std::string GetEccImpl::getConfigIds(const Ice::Current &)
{
	rc::ConfigIdSet configIdSet;
	try
	{

		backEnd().getConfigIds(configIdSet);
	}
	catch (const get::rc::SM::Exception & e)
	{
		throw mdaq::utl::CmdException(e.errorMessage);
	}
	return configIdSet.toXml();
}
//__________________________________________________________________________________________________
std::string GetEccImpl::getSubConfigIds(const std::string & subConfigType, const Ice::Current &)
{
	rc::SubConfigIdSet subConfigIds;
	backEnd().getSubConfigIds(subConfigType, subConfigIds);
	return subConfigIds.toXml();
}
//__________________________________________________________________________________________________
void GetEccImpl::loadPLLConfig(const Ice::Float ckwFreq_MHz, const Ice::Current &)
{
	backEnd().loadPLLConfig(ckwFreq_MHz);
}
//__________________________________________________________________________________________________
void GetEccImpl::writePLLRegister(const Ice::Long value, const Ice::Current &)
{
	backEnd().writePLLRegister(value);
}
//__________________________________________________________________________________________________
void GetEccImpl::resetPLL(const Ice::Current &)
{
	backEnd().resetPLL();
}
//__________________________________________________________________________________________________
/**
 * 	Sets up CoBo-Mutant link mode:
 *	  - Normal Operation (Input buffer enabled)
 *	  - CoBo Disconnected from Mutant network
 *	  - GMC Feedback Mode (Alignment)
 *	  - WSCA Feedback Mode (Alignment)
 */
void GetEccImpl::setupCoBoMutantLink(get::cobo::MutantLinkMode mode, const Ice::Current &)
{
	backEnd().setupCoBoMutantLink(mode);
}
//__________________________________________________________________________________________________
/**
 * Return reference to MuTanT alignment server.
 * @param c Information about the current method invocation.
 * @return Proxy to alignment server of current MuTanT node.
 */
mutant::AlignmentServerPrx GetEccImpl::mutantAlignmentServer(const Ice::Current & /* c */)
{
	return backEnd().mutantAlignmentServer();
}
//__________________________________________________________________________________________________
} // get
