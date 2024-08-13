/* =====================================================================================================================
 * EccClient.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class EccClient
 * Created on: Jan 5, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 *               Frederic Druillole (frederic.druillole@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GET (General Electronics for TPC) data acquisition framework @ CEA-Irfu.
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

#include "get/EccClient.h"
#include <CCfg/CConfig.h>
#include <CCfg/Io/BufferStorage.h>
#include <CCfg/Io/Document.h>
#include <CCfg/Io/IncludeHandler.h>
#include <CCfg/Io/XmlFormat.h>
#include "get/asad/utl.h"
#include "mdaq/utl/Server.h"
#include "utl/BitFieldHelper.hpp"
#include "utl/Logging.h"
#include <boost/lexical_cast.hpp>
#include <iomanip>
using mdaq::utl::CmdException;
#include <QFile>
#include <QFileInfo>
#include <QString>

#define DBG_IF(cond)  if(!(cond)) ; else LOG_DEBUG()
#define DEBUG_START()  DBG_IF(verbose) << "&gt;&gt;&gt;&gt; "
#define DEBUG_END()    DBG_IF(verbose) << "&lt;&lt;&lt;&lt; "
#undef MSG_IF
#define MSG_IF(cond)  DBG_IF(cond)

/**
 * Static variables
 */
bool EccClient::verbose = true;

const uint16_t EccClient::KNOWN_AGET_VERSIONS[NUM_KNOWN_AGET_VERSIONS] = {0x201, 0x202, 0x203}; /// ASTRE = 0x203

//__________________________________________________________________________________________________
EccClient::EccClient(const std::string& ipDotDecimalAndPort, const std::string& servantName)
	: sockAddr(ipDotDecimalAndPort), agetMask_(0), asadMask_(0), userCalibChipMask_(0), servantName(servantName)
{
}
//__________________________________________________________________________________________________
get::GetEccPrx EccClient::ecc()
{
	if (!ecc_)
	{
		const std::string proxyString = sockAddr.iceProxyString(servantName);
		LOG_INFO() << "Getting ECC proxy: " << proxyString;
		Ice::CommunicatorPtr commmunicator = ::mdaq::utl::Server::ic();
		try
		{
			ecc_ = get::GetEccPrx::checkedCast(commmunicator->stringToProxy(proxyString));
		}
		catch (const Ice::ConnectFailedException & e)
		{
			throw CmdException(std::string("Could not connect to ECC at '") + proxyString + "': " + e.what());
		}

		if (!ecc_)
		{
			throw CmdException(std::string("Could not connect to ECC at '") + proxyString + "'!");
		}
	}
	return ecc_;
}
//__________________________________________________________________________________________________
/**
 * Returns ID (0 to 15) of first active chip
 */
size_t EccClient::firstActiveChip() const
{
	size_t id=0;
	for (size_t id=0; id<16; ++id)
	{
		if (utl::BitFieldHelper< uint16_t >::getBit(agetMask_, id))
			break;
	}
	return id;
}
//__________________________________________________________________________________________________
/**
 * Connects to register server with given IP address and port.
 * @param targetAddress Target socket address (IP address and port).
 */
void EccClient::connectNode(const std::string & targetAddress)
{
	ecc()->connect("node", targetAddress);
}
//__________________________________________________________________________________________________
/**
 * Initializes hardware for node with given IP address according to configuration whose path is given.
 * @param targetAddress Target socket address.
 * @param hwPath Path of hardware configuration.
 */
void EccClient::loadHwDescription(const std::string & targetAddress, const std::string & hwPath)
{
	DEBUG_START() << "loadHwDescription(" << targetAddress << ", " << hwPath << ')';

	LOG_INFO() << "Connecting to node at " << targetAddress;
	ecc()->removeAllNodes();
	connectNode(targetAddress);

	LOG_INFO() << "Loading '" << QFileInfo(QString::fromStdString(hwPath)).fileName().toStdString() << '\'';
	CCfg::Io::Document doc;
	// Read hardware description
	try
	{
		doc.loadFromFile(hwPath);
	}
	catch (const CCfg::Io::Document::IoError & e)
	{
		throw mdaq::utl::CmdException(std::string("Error reading hardware description file.") + e.what());
	}

	// Store description into an XML string
	CCfg::Io::BufferStorage buffer;
	CCfg::Io::XmlFormat xml;
	xml.format(doc.getConfig(), buffer, CCfg::Io::NoOpIncludeHandler::getInstance());
	std::string hwDescr;
	buffer.readAll(hwDescr);

	// Load hardware description
	ecc()->describeNode("CoBo", "0", hwDescr);

	DEBUG_END() << "loadHwDescription";
}
//__________________________________________________________________________________________________
/**
 * Creates embedded software connection to GetController DAQ server.
 */
void EccClient::daqConnect(const std::string & dataRouterAddress, const std::string & flowType)
{
	DEBUG_START() << "daqConnect";
	try
	{
		DBG_IF(verbose) << "Connecting data sender to router of type '" << flowType << "' at " << dataRouterAddress;
		std::string objectName = std::string("daq-") + flowType;
		ecc()->connect(objectName, dataRouterAddress);
		std::cout << " Data sender successfully created.                  " << std::endl;
		LOG_INFO() << "Data sender successfully created.";
	}
	catch (const CmdException & e)
	{
		LOG_ERROR() << "Could not create data sender!";
		throw;
	}
	DEBUG_END() << "daqConnect";
}
//__________________________________________________________________________________________________
void EccClient::daqDisconnect()
{
	LOG_DEBUG() << "Disconnecting from data router...";
	ecc()->connect("daq-none", "");
}
//__________________________________________________________________________________________________
/**
 * Checks hardware version.
 *
 * The CoBo register 'hardwareVersion' (offset 0x1D4) was introduced with firmware 20131107.
 * It has two valid values: 1= Reduced CoBo, 2= CoBo >= v0.2.
 *
 * @return Value of hardware version register or 1 if not found.
 */
uint32_t EccClient::checkHardwareVersion()
{
	ecc()->connect("device", "ctrl");

	std::ostringstream oss;
	oss << "Checking hardware version... ";
	uint32_t hardwareVersion = 1; // Assume R-CoBo if register is not described
	try
	{
		hardwareVersion = ecc()->readRegister("hardwareVersion");
		std::string hardwareVersionStr = "R-CoBo";
		if (2 <= hardwareVersion)
		{
			hardwareVersionStr = "CoBo >= v0.2";
		}
		else
		{
			hardwareVersion = 1;
		}
		oss << hardwareVersion << " (" << hardwareVersionStr << ')';
	}
	catch (const CmdException &)
	{
		oss << " not found";
	}
	LOG_INFO() << oss.str();
	return hardwareVersion;
}
//__________________________________________________________________________________________________
/**
 * Checks firmware release  number and release date.
 *
 * They are provided by CoBo registers 'firmwareDate' (offset 0x1D8) and 'firmwareRelease' (offset 0x1DC).
 * firmwareDate is a hexadecimal value displaying the build date of the firmware in the format 0xYYYYMMDD (ie. firmware from Nov 7th, 2013 will be value 0x20131107).
 * firmwareRelease stores the sequential release number of the firmware.  The firmware build on Nov 7th, 2013 will be version 1 since it is the first version with this register.
 */
uint32_t EccClient::checkFirmwareRelease()
{
	ecc()->connect("device", "ctrl");

	std::ostringstream oss;
	oss << "Checking firmware release... ";
	uint32_t firmwareRelease = 0;
	try
	{
		firmwareRelease = ecc()->readRegister("firmwareRelease");
		const uint32_t firmwareDate = ecc()->readRegister("firmwareDate");
		oss << firmwareRelease << " (" << std::hex << firmwareDate << ')';
	}
	catch (const CmdException &)
	{
		oss << " not found";
	}
	LOG_INFO() << oss.str();
	return firmwareRelease;
}
//__________________________________________________________________________________________________
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
void EccClient::setCoBoInitDone(bool done)
{
	try
	{
		ecc()->connect("device", "ctrl");
		ecc()->writeField("mutantConfig", "initDone", done);
	}
	catch (const Ice::LocalException & e)
	{
		LOG_WARN() << "Could not set CoBo 'initDone' bit to " << done << " (" << e.ice_name() << ')';
	}
	catch (const mdaq::utl::CmdException & e)
	{
		LOG_WARN() << "Could not set CoBo 'initDone' bit to " << done << ". " << e.reason << ')';
	}
}
//__________________________________________________________________________________________________
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
void EccClient::setCoBoMutantMode(bool enabled)
{
	try
	{
		ecc()->connect("device", "ctrl");
		ecc()->writeField("mutantConfig", "mutantMode", enabled);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		LOG_WARN() << e.what();
	}
}
//__________________________________________________________________________________________________
/**
 * Sets the configuration bits of the data flow
 * to their initial values. See script beforeStart.ecc.
 */
void EccClient::pipeInit()
{
	LOG_INFO() << "Initializing CoBo firmware pipeline";
	DEBUG_START() << "pipeInit";

	ecc()->connect("device", "ctrl");

	/* bits of the pipeCtrl register */
	ecc()->writeField("pipeCtrl", "multiplicityReadOut", 0);	// Disable creation of frames with multiplicity data (depreciated, specific to 1 AGET).
	ecc()->writeField("pipeCtrl", "modeEnable_2p", 0);			// Disable 2p mode.
	ecc()->writeField("pipeCtrl", "bitSlipCalibrate", 0);		// Disable CoBo deserializer calibration sequence.
	ecc()->writeField("pipeCtrl", "partialReadOut", 0);			// Disable partial readout mode.
	ecc()->writeField("pipeCtrl", "testMode", 0);				// Disable test mode (where frames contain incremental counter instead of actual data).
	ecc()->writeField("pipeCtrl", "zeroSuppressionEnable", 0);	// Disable zero suppression mode.
	ecc()->writeField("pipeCtrl", "firmwareEnable", 0);			// Disable the firmware state machine.
	ecc()->writeField("pipeCtrl", "hitRegisterOverride", 0);	// Depreciated.

	// Corresponds to script sequencer.ecc
	ecc()->writeField("pipeCtrl", "SCW", 0);
	ecc()->writeField("pipeCtrl", "SCR", 0);
	ecc()->writeField("pipeCtrl", "ScrScwOverride", 1);			// Allow software override of the SCW and SCR bit fields.

	ecc()->writeField("pipeCtrl", "reset", 1);					// Reset all the CoBo logic.
	ecc()->writeField("pipeCtrl", "reset", 0);
	ecc()->writeField("pipeCtrl", "backpressure", 0);			// Depreciated (specific to 1 AGET).
	ecc()->writeField("pipeCtrl", "intEnable", 0);
	DEBUG_END() << "pipeInit";
}
//__________________________________________________________________________________________________
/**
 * Sets the configuration bits of the trigger core to their initial values.
 */
void EccClient::triggerInit()
{
	LOG_INFO() << "Initializing trigger core";
	DEBUG_START() << "triggerInit";

	ecc()->connect("device", "ctrl");

	// Set manual delay between rising edge of SCR and first data on CoBo
	// 15 corresponds to a 1 meter cable
	// 16 corresponds to a 3 meter cable (not sure)
	for (size_t agetId=0; agetId < 4u; agetId++)
	{
		std::string fieldName = std::string("AGET") + boost::lexical_cast< std::string >(agetId);
		ecc()->writeField("readDataDelay", fieldName, 14);
	}

	// Delay (units of 10 ns) between event causing the trigger and actual freeze of the AGET memory.
	ecc()->writeRegister("triggerDelay", 0);

	// Delay (units of 10 ns) between event causing the second trigger and actual freeze of the AGET second memory.
	ecc()->writeRegister("triggerDelay_2p", 0);

	// Delay (units of 10 ns) after the first trigger after which the second trigger can be issued.
	ecc()->writeRegister("triggerWindowStart_2p", 0);

	// Maximum duration (units of 10 ns) to wait for the second trigger. After this delay, first trigger is canceled.
	ecc()->writeRegister("triggerWindowEnd_2p", 100000); // 1 ms

	// Multiplicity threshold for the second trigger.
	ecc()->writeRegister("multThreshold_2p", 3000);

	// Period of the second trigger if periodic, in units of 10 ns.
	ecc()->writeRegister("triggerTime_10ns_2p", 1000); // 10 us

	// Multiplicity threshold (total of the 4 AsAd's).
	ecc()->writeRegister("multThreshold", 3000);

	// Period of the trigger, if periodic, in units of 10 ns.
	ecc()->writeRegister("triggerPeriod_10ns", 500000000); // 5s

	// Trigger mode
	ecc()->writeRegister("triggerMask", 0); // no trigger mode selected

	// Number of samples to sum in the multiplicity sliding window for the multiplicity mode.
	ecc()->writeRegister("multWindSize", 16);

	// Number of time buckets to readout.
	ecc()->writeRegister("readDepth", 512); // is set by default to the reading of all time buckets

	// Channels to always read in partial readout mode
	ecc()->writeRegister("readAlwaysMask1_4", 0);
	ecc()->writeRegister("readAlwaysMask5_36", 0);
	ecc()->writeRegister("readAlwaysMask37_68", 0);
	// Channels to read only if hit in partial readout mode.
	ecc()->writeRegister("readIfHitMask1_4", 0xF);
	ecc()->writeRegister("readIfHitMask5_36", 0xFFFFFFFF);
	ecc()->writeRegister("readIfHitMask37_68", 0xFFFFFFFF);

	DEBUG_END() << "triggerInit";
}
//__________________________________________________________________________________________________
/**
 * Sets the configuration bits of the circular buffer
 * to their initial values
 */
void EccClient::circularBufferInit()
{
	LOG_INFO() << "Initializing circular buffers";
	DEBUG_START() << "circularBufferInit";

	ecc()->connect("device", "ctrl");

	// Memory addresses of the circular buffers of the 4 AsAd boards
	// Total space = 0x8000000 to 0x10000000 = 128 MiB ?
	// FIXME: If only 1 AsAd board is used, it should be able to use all 128 MiB
	uint64_t memStart_asad0 = UINT64_C(0x8000000);
	uint64_t bufferSize_B = UINT64_C(0x2000000); // 32 MiB
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		uint64_t memStart = memStart_asad0 + asadIdx*bufferSize_B;
		// Start and end
		ecc()->writeRegister(std::string("memStart_asad") + boost::lexical_cast< std::string >(asadIdx), memStart);
		ecc()->writeRegister(std::string("memEnd_asad") + boost::lexical_cast< std::string >(asadIdx), memStart + bufferSize_B);
		// Read pointer
		ecc()->writeRegister(std::string("readPtr_asad") + boost::lexical_cast< std::string >(asadIdx), memStart);
	}

	// Number of events per CPU interrupt
	try
	{
		ecc()->writeRegister("interruptRate", 0x1);
	}
	catch (const CmdException &) {}

	// Threshold (Bytes) indicating the amount of free space below which the firmware stops writing in the circular buffer.
	// Should be greater than the maximum frame size.
	ecc()->writeRegister("aFullLevel", 0x01000000);

	// Level (Bytes) below which the firmware resumes writing data to the circular buffer (see 'aFullLevel').
	ecc()->writeRegister("aEmptyLevel", 0x01000000);

	DEBUG_END() << "circularBufferInit";
}
//__________________________________________________________________________________________________
/**
 * Sets the multiplicity readout mode to get multiplicity analog value.
 * Specific to 1 AGET firmware. Creates frames of multiplicity data.
 */
void EccClient::coboMultReadMode(bool enable)
{
	DEBUG_START() << "coboMultReadMode";
	ecc()->connect("device", "ctrl");
	if (enable) ecc()->writeField("pipeCtrl", "firmwareEnable", 1);
	ecc()->writeField("pipeCtrl","multiplicityReadOut",enable); // 100 MHz
	if (not enable) ecc()->writeField("pipeCtrl", "firmwareEnable", 0);
	DEBUG_END() << "coboMultReadMode";
}
//__________________________________________________________________________________________________
/**
 * Checks whether AsAd board is active.
 */
bool EccClient::isAsAdBoardActive(const size_t & asadIdx) const
{
	return ::utl::BitFieldHelper< uint16_t >::getBit(asadMask_, asadIdx);
}
/**
 * Read AGET version register.
 */
uint16_t EccClient::agetVersion(size_t asadId, size_t agetId)
{
	// Select AGET to read from
	selectAgetForRead(asadId, agetId);
	// Read AGET register 5
	ecc()->connect("device", "aget");
	return ecc()->readRegister("reg5");
}
//__________________________________________________________________________________________________
void EccClient::clearAllLeds()
{
	try
	{
		ecc()->ledManager()->setLEDs(false);
	}
	catch (const Ice::OperationNotExistException & e)
	{
		LOG_WARN() << e.what();
	}
}
//__________________________________________________________________________________________________
void EccClient::setLedState(const LedType & type, const size_t & asadId, const LedState & state)
{
	try
	{
		ecc()->ledManager()->modifyLED(type, asadId, state);
	}
	catch (const Ice::OperationNotExistException & e)
	{
		LOG_WARN() << e.what();
	}
}
//__________________________________________________________________________________________________
/**
 * Powers off all four AsAd boards.
 */
void EccClient::asadPowerOff()
{
	LOG_INFO() << "Powering off AsAd board(s)";
	try
	{
		// Check monitoring status of active AsAd boards (it will no longer be possible after they are powered off).
		ecc()->connect("device", "ctrl");
		uint32_t asadEnable = ecc()->readRegister("asadEnable");
		for (size_t asadId = 0; asadId < 4u; asadId++)
		{
			if (not ::utl::BitFieldHelper< uint32_t >::getBit(asadEnable, asadId)) continue;
			selectAsAdForRead(asadId);
			const uint16_t status = getAsAdAlertStatus();
			if (0 != status)
			{
				LOG_WARN() << "AsAd board no. " << asadId << " is in auto-protection mode (status "
				        << std::hex << std::showbase << status << std::dec << " e.g. "
				        << get::asad::utl::stringifyMonitoringAlertStatus(status) << ')';
			}
		}

		// Power off all AsAd boards
		ecc()->connect("device", "ctrl");
		ecc()->writeRegister("asadEnable", 0);
	}
	catch (::mdaq::utl::CmdException &)
	{
	}
}
//__________________________________________________________________________________________________
/**
 * Powers ON AsAd board(s) and sets the configuration bits
 * of the AsAd slow controller and of AsAd to their initial values.
 */
void EccClient::asadInit(int time_ms, const uint8_t & abortIfNoPowerSupply)
{
	LOG_INFO() << "Initializing AsAd boards";

	clearAllLeds();

	// Resets the AsAd slow controller on CoBo
	ecc()->connect("device", "ctrl");
	ecc()->writeField("asadController", "reset", 1); // If pulsed, resets error flag, busy flag and checksum error
	ecc()->writeField("asadController", "manualCal", 0); // Value of the manual delay (units of 10 ns)
	ecc()->writeField("asadController", "calibration", 0); // If set, calibration is automatic (1 bit per AsAd)
	ecc()->writeField("asadController", "reset", 0);

	// Powers ON AsAd boards
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		// Acknowledge AsAd alert status
		// As long as the alert status register content is not readout the VDD voltage regulator cannot be switched on again.
		if (isAsAdConnected(asadIdx))
		{
			const uint16_t mask = (1ul << asadIdx);
			setAsAdSlowControlMask(mask);
			selectAsAdForRead(asadIdx);
			// Read alert status
			const uint16_t alertStatus = getAsAdAlertStatus();
			if (0 != alertStatus)
			{
				LOG_WARN() << "AsAd board no. " << asadIdx << " was in auto-protection mode (status "
						<< std::hex << std::showbase << alertStatus << std::dec
						<< " e.g. " << get::asad::utl::stringifyMonitoringAlertStatus(alertStatus) << ')';
			}
		}

		// Power OFF board
		std::string fieldName = std::string("powerON") + boost::lexical_cast< std::string >(asadIdx);
		ecc()->connect("device", "ctrl");
		ecc()->writeField("asadEnable", fieldName, 0);
		sleep_ms(time_ms);
		if (isAsAdBoardActive(asadIdx))
		{
			// Check board is connected
			if (not isAsAdConnected(asadIdx))
			{
				setLedState(get::cobo::LedB, asadIdx, get::cobo::LedSlowPulse);
				throw CmdException(std::string("AsAd board no. ") + boost::lexical_cast< std::string >(asadIdx) + std::string(" is not connected to CoBo!"));
			}
			// Check board is power supplied
			if (not isAsAdPowerSupplied(asadIdx))
			{
				std::string errorMsg = std::string("AsAd board no. ") + boost::lexical_cast< std::string >(asadIdx) + std::string(" seems to be disconnected from its power supply!");
				if (::utl::BitFieldHelper< uint8_t >::getBit(abortIfNoPowerSupply, asadIdx))
				{
					setLedState(get::cobo::LedB, asadIdx, get::cobo::LedSlowPulse);
					throw CmdException(errorMsg);
				}
				else
					LOG_ERROR() << errorMsg;
			}
			// Power ON board
			DBG_IF(verbose) << "Powering on AsAd board no. " << asadIdx;
			ecc()->writeField("asadEnable", fieldName, 1);
			sleep_ms(time_ms);

			// Check board status
			bool plg, alarm, pwu;
			asadPowerStatus(asadIdx, plg, alarm, pwu);
			if (alarm)
			{
				std::ostringstream msg;
				msg << "AsAd board no. " << asadIdx << " was not successfully powered on"
						" (PLG=" << (not plg) << ", AL=" << ", PWU=" << pwu << "!";
				setLedState(get::cobo::LedB, asadIdx, get::cobo::LedSlowPulse);
				throw mdaq::utl::CmdException(msg.str());
			}

			// Switch 'P' LED on to indicate AsAd board is successfully connected and powered.
			setLedState(get::cobo::LedP, asadIdx, get::cobo::LedOn);
		}
	}
}
//__________________________________________________________________________________________________
/**
 * Initializes AsAd DAC device registers.
 * It seems the value of the DAC1 ("DATA") register is not always its official power-on value of 0.
 * @see GET-AS-002-0010.
 */
void EccClient::asadDacInit()
{
	ecc()->writeRegister("DAC0", 0x81);
	ecc()->writeRegister("DAC1", 0);
	// Register DAC2 (VERSION) is read only.
	ecc()->writeRegister("DAC3", 0);
	ecc()->writeRegister("DAC4", 0);
	ecc()->writeRegister("DAC5", 0);
	ecc()->writeRegister("DAC6", 0); // 0x38 according to GET-AS-002-0010
	ecc()->writeRegister("DAC7", 0);
}
//_________________________________________________________________________________________________
void EccClient::asadPowerStatus(const size_t & asadIdx, bool & plg, bool & al, bool & pwu)
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
uint32_t EccClient::asadGetIdentifier()
{
	uint32_t identifier = 0;
	try
	{
		// Read AIM registers one by one, in the right order
		ecc()->connect("device", "asad");
		for (size_t i=0; i < 4; ++i)
		{
			const std::string regName = std::string("AIM") + boost::lexical_cast< std::string >(i);
			const uint8_t regValue = ecc()->readRegister(regName);
			utl::BitFieldHelper< uint32_t >::setField(identifier, 8u*(3u-i), 8u, regValue);
		}
	}
	catch (CmdException & e)
	{
		;
	}
	return identifier;
}
//__________________________________________________________________________________________________
/**
 * Sets AsAd board user ID.
 *
 * The AsAd identifier is a four 8-bit register (or a 32-bit register, if you prefer) composed of
 * - three 8-bit "read -only" fields (serial number etc...)
 * - and one 8-bit "read-write" field that allows the user to define its own label to identify a particular AsAd.
 */
void EccClient::asadSetBoardUserID(const uint8_t & userID)
{
	ecc()->connect("device", "asad");
	ecc()->writeRegister("AIM0", userID);
}
//__________________________________________________________________________________________________
/**
 * Initializes 2p mode readout.
 */
void EccClient::aget2pInitNow()
{
	MSG_IF(verbose) << "Initializing 2p mode readout";
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("aget_2p_initNow", 1);
	sleep_ms(20);
	ecc()->writeRegister("aget_2p_initNow", 0);
}
//__________________________________________________________________________________________________
/**
 * Sets mask of active AGET chips.
 */
void EccClient::setActiveAGetMask(const uint16_t mask)
{
	agetMask_ = mask;
	LOG_INFO() << "Mask of active AGET chips: " << std::showbase << std::hex << agetMask_ << std::dec;
}
//__________________________________________________________________________________________________
/**
 * Sets mask of active AsAd boards.
 */
void EccClient::setActiveAsAdMask(const uint16_t mask)
{
	asadMask_ = mask;
	LOG_INFO() << "Mask of active AsAd boards: " << std::showbase << std::hex << asadMask_ << std::dec;
}
//__________________________________________________________________________________________________
void EccClient::selectCalibChips(const uint8_t calibChipMask)
{
	userCalibChipMask_ = calibChipMask;
}
//__________________________________________________________________________________________________
/**
 * Initializes the AGET chip, together with its slow control
 */
void EccClient::agetInit()
{
	DEBUG_START() << "agetInit";

	ecc()->connect("device", "ctrl");

	// Reset AGET controller error and busy flags
	ecc()->writeField("agetController", "reset", 0);
	ecc()->writeField("agetController", "reset", 1);
	ecc()->writeField("agetController", "reset", 0);

	// Set mask for broadcast of AGET slow control commands
	setAgetSlowControlMask(agetMask_);
	ecc()->readField("aget_sc_configuration", "SC_En");

	// Set mask for broadcast of AsAd slow control commands
	setAsAdSlowControlMask(asadMask_);

	// Reset the AsAd state machine for the AGET slow control
	ecc()->connect("device", "asad");
	ecc()->writeRegister("resetFast", 0);  // Reset internal state machine of Fast Control process

	// Calibrate AGET slow control
	agetCal();

	ecc()->connect("device", "aget");
	ecc()->writeField("reg1", "DAC", 1);
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("initSCAReadOutAfterPowerOnNow", 1);

	DEBUG_END() << "agetInit";
}
//__________________________________________________________________________________________________
/**
 * Triggers the automatic calibration of the data delay
 *  (time from rising edge of SCR signal until data appears at the CoBo).
 *
 * The calibration must be done after the other cable delays are calibrated
 *  and the local clocks manager of the AsAd has been configured.
 * It should be done before any registers of the AGET are configured,
 *  since AGET reg1 and reg2 are modified as part of the calibration.
 *
 * You can choose between the automatic or manual value by setting the SCRAutoDelayEnable fields in the pipeCtrl register.
 */
void EccClient::calibrateAutoReadDataDelay()
{
	LOG_INFO() << "Calibrating SCR data delays";
	// Based on script dataCal.ecc
	ecc()->connect("device", "ctrl");
	ecc()->writeField("aget_sc_configuration", "SC_En", agetMask_);

	ecc()->connect("device", "aget");
	ecc()->writeRegister("reg1", 0x3e3a41a7);
	ecc()->writeRegister("reg2", 0x0028a197);

	ecc()->connect("device", "ctrl");
	try
	{
		ecc()->writeField("pipeCtrl", "SCRAutoDelayEnable", 1);
	}
	catch (const CmdException &)
	{
		ecc()->writeField("pipeCtrl", "SCRAutoDelayEnable0", 1); // Field renamed on April 11th, 2013
	}
	ecc()->writeField("pipeCtrl", "SCRAutoDelayCalibrate", 1);
	ecc()->writeField("pipeCtrl", "intEnable", 0);
	ecc()->writeField("pipeCtrl", "SCR", 0);
	ecc()->writeField("pipeCtrl", "SCW", 0);
	ecc()->writeField("pipeCtrl", "ScrScwOverride", 1);
	ecc()->writeField("pipeCtrl", "reset", 0);
	ecc()->writeField("pipeCtrl", "firmwareEnable", 1);
	ecc()->writeField("pipeCtrl", "SCW", 1);
	sleep_ms(100);
	ecc()->writeField("pipeCtrl", "SCW", 0);
	ecc()->writeField("pipeCtrl", "SCR", 1);
	ecc()->writeField("pipeCtrl", "SCR", 0);
	sleep_ms(100);
	ecc()->writeField("pipeCtrl", "firmwareEnable", 0);
	ecc()->writeField("pipeCtrl", "SCRAutoDelayCalibrate", 0);
	ecc()->writeField("pipeCtrl", "ScrScwOverride", 0);
}
//__________________________________________________________________________________________________
/**
 * Enables or disables the use of automatically calibrated read data delays.
 *
 * Sets the SCRAutoDelayEnable CoBo register.
 * Ignored in the microTCA CoBo firmware?
 * See email from Nathan Usher on April 11th and August 2nd 2013.
 */
void EccClient::asadSetAutoReadDataDelayEnabled(const size_t & asadId, bool enabled)
{
	LOG_INFO() << "Selecting " << (enabled?"automatic":"manual") << " data delays for AsAd no. " << asadId;
	ecc()->connect("device", "ctrl");
	std::string fieldName = std::string("SCRAutoDelayEnable");
	try
	{
		ecc()->writeField("pipeCtrl", fieldName, enabled ? 1 : 0);
	}
	catch (const CmdException &)
	{
		fieldName = std::string("SCRAutoDelayEnable") + boost::lexical_cast<std::string>(asadId);
		ecc()->writeField("pipeCtrl", fieldName, enabled ? 1 : 0); // Field renamed on April 11th, 2013
	}
}
//_________________________________________________________________________________________________
/**
 * Set AGET input protections.
 * Two goals:
 * 	- Enable AsAd ZAP protections since AsAs is now totally powered.
 * 	- Prevent AGET sampling or readout commands until clocks are prepared.
 *
 */
void EccClient::asadSetAgetInputProtections()
{
	LOG_DEBUG() << "Setting AGET input protections";

	// Write to all AsAd boards simultaneously
	setAsAdSlowControlMask(asadMask_);

	ecc()->connect("device", "asad");
	// Configure connection to TPC pads and anti-sparkling protection
	ecc()->writeRegister("IM0", 0x0F);
	// Disable sampling and readout commands on some AGET inputs
	ecc()->writeRegister("IM1", 0x00);
}
//_________________________________________________________________________________________________
/**
 * Applies given configuration for AGET positions in AsAd.
 * Sets first register of the AsAd ADC device and the two registers of the AsAd Inputs Manager.
 */
void EccClient::asadConfigureAgetInputs(const uint8_t & adc0, const uint8_t & IM0, const uint8_t & IM1)
{
	ecc()->connect("device", "asad");
	DBG_IF(verbose) << "Setting ADC0 to " << std::showbase << std::hex << (uint16_t) adc0
			<< ", IM0 to " << (uint16_t) IM0 << " and IM1 to " << (uint16_t) IM1 << std::dec;
	// Power down some AGET's ADC channels
	ecc()->writeRegister("ADC0", adc0);
	// Configure connection to TPC pads and anti-sparkling protection
	ecc()->writeRegister("IM0", IM0);
	// Disable sampling and readout commands on some AGET inputs
	ecc()->writeRegister("IM1", IM1);
}
//_________________________________________________________________________________________________
/**
 * Runs AsAd ADC calibration sequence.
 * See email from N. Usher on June 13th, 2013 and Texas Instruments ADS6422 Data Sheet.
 */
void EccClient::asadAdcInit(ModeADC modeADC)
{
	LOG_INFO() << "Running AsAd ADC calibration sequence (" << (modeADC == DDR?"DDR":"SDR") << ')';
	// Set mask for broadcast of AsAd slow control commands
	setAsAdSlowControlMask(asadMask_);

	// Reset ADC
	adcReset();

	// Calibrate of the ADC lines
	adcCal(modeADC);

	// Put the ADC in its default configuration
	ecc()->connect("device", "asad");
	// Powers down ADC channels B, C, D
	//ecc()->writeRegister("ADC0", 0x1C);
	// Straight binary data format, normal ADC operation
	ecc()->writeRegister("ADC2", 0x200);
	// 2 Wire interface, bit clock selection
	ecc()->writeRegister("ADC5", (modeADC == SDR)?0x403:0x401);
}
//__________________________________________________________________________________________________
/**
 * Selects AsAd board(s) for slow control.
 */
void EccClient::setAgetSlowControlMask(const uint16_t & mask)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeField("aget_sc_configuration", "SC_En", mask);
}
//__________________________________________________________________________________________________
/**
 * Selects AsAd board(s) for slow control.
 */
void EccClient::setAsAdSlowControlMask(const uint16_t & mask)
{
	ecc()->connect("device", "ctrl");
	try
	{
		ecc()->writeField("ASAD_sc_mask", "select", mask);
	}
	catch (const CmdException & e)
	{
		LOG_WARN() << e.reason;
		LOG_WARN() << "You need to update your hardware description file!";
	}
}
//__________________________________________________________________________________________________
/**
 * Selects AGET chip to read from/
 */
void EccClient::selectAgetForRead(const size_t & asadId, const size_t & agetId)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeField("agetChipSelect", "d", 4*asadId + agetId);
}
//__________________________________________________________________________________________________
/**
 * Selects AsAd board(s) for read.
 */
void EccClient::selectAsAdForRead(const size_t & asadId)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeField("asadSelect", "d", asadId);
}
//__________________________________________________________________________________________________
/**
 * Selects a unique AsAd board for a slow control command.
 * @param asadIdx Index of the AsAd board to select.
 */
void EccClient::selectAsAdForSlowControl(const size_t & asadIdx)
{
	const uint16_t mask = (1ul << asadIdx);
	setAsAdSlowControlMask(mask);
	selectAsAdForRead(asadIdx);
}
//__________________________________________________________________________________________________
/**
 * Calibrates the AsAd slow control.
 */
void EccClient::asadCal()
{
	LOG_INFO() << "Calibrating AsAd slow control ...";

	setAsAdSlowControlMask(asadMask_);

	// Synchronizes the slow control lines
	ecc()->connect("device", "ctrl");
	ecc()->writeField("asadController", "reset", 1); // Reset error flag, busy flag and checksum error.
	ecc()->writeField("asadController", "reset", 0);
	ecc()->writeField("asadController", "manualCal", 0x0); // Value of the AsAd manual delay. Specific to Reduced CoBo.
	ecc()->writeField("asadController", "calibration", asadMask_); // If set, calibration is automatic. 1 bit per AsAd
	ecc()->writeRegister("asad_sc_start", 0x1);
	ecc()->writeField("asadController", "calibration", 0);

	// Check calibration of AsAd slow control was successful:
	// Manufacturer ID of monitoring device must be 0x41
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		// Skip unused board
		if (not isAsAdBoardActive(asadIdx)) continue;
		// Select AsAd board
		selectAsAdForSlowControl(asadIdx);
		// Read register
		ecc()->connect("device", "asad");
		uint64_t manufacturerID = ecc()->readRegister("monitorID");

		// Check manufacturer ID value is correct
		if (manufacturerID != UINT64_C(0x41))
		{
			std::ostringstream msg;
			msg << "Calibration of slow control of AsAd no. " << asadIdx << " failed: expected 0x41, read "
					<< std::hex << std::showbase << manufacturerID << '!';
			throw mdaq::utl::CmdException(msg.str());
		}
		// Read board identifier
		else
		{
			// Set user chosen bits of board identifier
			const uint8_t expectedUserID = 1 + asadIdx;
			asadSetBoardUserID(expectedUserID);

			// Read board identifier
			const uint32_t fullID = asadGetIdentifier();
			uint16_t actualUserID = utl::BitFieldHelper< uint32_t >::getField(fullID, 24u, 8u);
			uint32_t boardID = fullID & 0xFFFFFF;
			LOG_DEBUG() << "AsAd board no. " << asadIdx << " has identifier " << std::hex << boardID << std::dec
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
	LOG_INFO() << " ... AsAd slow control calibration done";
}
//__________________________________________________________________________________________________
/**
 * Triggers AsAd ADC software reset - resets all internal registers.
 */
void EccClient::adcReset()
{
	ecc()->connect("device", "asad");
	ecc()->writeRegister("ADC0", 0x400);
}
//__________________________________________________________________________________________________
/**
 * Sets which AGET should be used for cable delay calibrations.
 *
 * The 'calibrationAget' CoBo register (offset 0x1E0) sets which AGET should be used for cable delay calibrations.
 * It has 4 fields, one for each AsAd.  Each field is two bits.
 * This register must be set before the calibration sequence is started.
 * This setting should only be modified on AsAds that do not have an AGET installed in position 0.
 */
void EccClient::coboSetCalibrationAgetChips(const uint16_t & agetMask)
{
	LOG_INFO() << "Selecting calibration AGET chips (mask=" << std::hex << std::showbase << agetMask << std::dec << ')';
	try
	{
		// Loop over AsAd boards
		for (size_t asadIdx=0; asadIdx < 4; ++asadIdx)
		{
			size_t calibrationChipId = 0;
			if (isAsAdBoardActive(asadIdx))
			{
				// Find first active chip of the 4 chips on this AsAd board
				for (calibrationChipId=0; calibrationChipId < 4; ++calibrationChipId)
				{
					if (utl::BitFieldHelper< uint16_t >::getBit(agetMask, asadIdx*4+calibrationChipId))
						break;
				}
				// Check which chip the user wants to use
				size_t userCalibrationChipId = utl::BitFieldHelper< uint8_t >::getField(userCalibChipMask_, asadIdx*2, 2);
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
				LOG_INFO() << "Selecting AGET chip no. " << asadIdx*4+calibrationChipId << " for AsAd board no. " << asadIdx << " serial control delay calibration";
			}
			// Set AGET chip for serial control delay calibration
			ecc()->connect("device", "ctrl");
			std::string fieldName = std::string("AsAd") + boost::lexical_cast< std::string >(asadIdx);
			//LOG_DEBUG() << "Setting calibrationAget." << fieldName << " to " << calibrationChipId;
			ecc()->writeField("calibrationAget", fieldName, calibrationChipId);
		}
	}
	catch (const CmdException & e)
	{
		// Case of old firmware/hardware description where the 'calibrationAget' register does not exist.
		LOG_WARN() << e.what();
	}
}
//__________________________________________________________________________________________________
/**
 * Calibrates the AGET slow control.
 *
 * By default, during the calibration sequence, the firmware measures the time delay for data received from AGET 0.
 */
void EccClient::agetCal()
{
	LOG_INFO() << "Calibrating AGET slow control";
	DEBUG_START() << "agetCal";

	// Create local mask for easier debugging
	const uint16_t mask = agetMask_;

	// Select which AGET to use for serial control delay calibration
	coboSetCalibrationAgetChips(agetMask_);

	adcReset();
	adcTest(DDR);

	ecc()->connect("device", "ctrl");
	ecc()->writeField("ckpFastDivider", "d", 0);
	ecc()->writeField("agetController", "reset", 1);
	ecc()->writeField("agetController", "reset", 0);

	ecc()->writeField("aget_sc_configuration", "SC_En", mask);
	ecc()->writeField("aget_sc_configuration", "pulseEnable", mask);

	ecc()->writeField("agetController", "hitRegister", 0);
	ecc()->writeField("agetController", "calibration", 0);
	//ecc()->writeRegister("aget_sc_rw", 1); // < from agetCal.ecc

	ecc()->connect("device", "aget");
	ecc()->readRegister("reg5");
	ecc()->connect("device", "ctrl");
	ecc()->writeField("agetController", "reset", 1);
	ecc()->writeField("agetController", "reset", 0);

	ecc()->writeField("aget_sc_configuration", "pulseEnable", mask);

	ecc()->writeField("agetController", "hitRegister", 0);
	ecc()->writeField("agetController", "calibration",1);

	ecc()->connect("device", "aget");
	ecc()->readRegister("reg5");
	ecc()->connect("device", "ctrl");
	ecc()->writeField("agetController", "reset", 1);
	ecc()->writeField("agetController", "reset", 0);
	ecc()->writeField("agetController", "calibration", 0);

	ecc()->writeField("aget_sc_configuration", "pulseEnable", mask);

	ecc()->writeField("agetController", "hitRegister", 1);
	ecc()->connect("device", "aget");
	ecc()->writeRegister("hit", UINT64_C(0xA000000000000005));
	ecc()->connect("device", "ctrl");
	ecc()->writeField("agetController", "calibration", 1);

	ecc()->writeField("aget_sc_configuration", "pulseEnable", 0);

	ecc()->writeField("agetController", "hitRegister", 1);
	ecc()->connect("device", "aget");
	ecc()->writeRegister("hit", UINT64_C(0xA000000000000005));
	ecc()->connect("device", "ctrl");
	ecc()->writeField("agetController", "calibration", 0);

	ecc()->writeField("aget_sc_configuration", "pulseEnable", mask);
	ecc()->writeField("agetController", "hitRegister", 0);


	ecc()->connect("device", "aget");
	ecc()->readRegister("reg5");

	ecc()->connect("device", "ctrl");
	ecc()->writeField("agetController", "reset", 1);
	ecc()->writeField("agetController", "reset", 0);

	// Print AGET chip version numbers so that user can see if calibration of AEGT slow control was successful
	printAgetChipVersionNumbers();

	DEBUG_END() << "agetCal";
}
//__________________________________________________________________________________________________
/**
 * Print version numbers of active AGET chips
 */
void EccClient::printAgetChipVersionNumbers(bool check)
{
	// Switch to slow control mode to be able to read AGET register 5
	switchToSlowControl();

	if (check)
	{
		LOG_DEBUG() << "Checking AGET chip version numbers are valid...";
	}

	std::ostringstream msg;
	msg << "AGET chip version numbers:";
	for (size_t asadId = 0; asadId < 4u; ++asadId)
	{
		for (size_t agetId=0; agetId < 4u; ++agetId)
		{
			if (utl::BitFieldHelper< uint16_t >::getBit(agetMask_, 4*asadId+agetId))
			{
				selectAgetForRead(asadId, agetId);
				ecc()->connect("device", "aget");
				const uint16_t agetVersion = (uint16_t) ecc()->readRegister("reg5");
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
	LOG_DEBUG() << msg.str();
}
//__________________________________________________________________________________________________
/**
 * Calibration of the AsAd ADC.
 * See email from N. Usher on July 13th, 2013 and Texas Instruments ADS6422 Data Sheet.
 */
void EccClient::adcCal(ModeADC modeADC)
{
	ecc()->connect("device", "asad");
	// Straight binary data format, DESKEW test pattern (serial stream of 1010..)
	ecc()->writeRegister("ADC2", 0x2C0);
	// Custom test pattern = 0xC3 ?
	ecc()->writeRegister("ADC3", 0x0C3);
	ecc()->writeRegister("ADC4", 0x000);
	// 2 Wire interface, bit clock selection
	ecc()->writeRegister("ADC5", (modeADC == SDR)?0x403:0x401);

	// Initiate the CoBo deserializer deskew sequence
	ecc()->connect("device", "ctrl");
	try
	{
		ecc()->writeField("pipeCtrl", "adcDeskew",1);
		IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(1000));
		ecc()->writeField("pipeCtrl", "adcDeskew", 0);
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		LOG_WARN() << "You need to update your hardware description file!";
	}

	ecc()->connect("device", "asad");
	// Straight binary data format, CUSTOM test pattern (registers ADC3 and ADC4)
	ecc()->writeRegister("ADC2", 0x2A0);

	// Initiate the CoBo deserializer calibration sequence
	ecc()->connect("device", "ctrl");
	ecc()->writeField("pipeCtrl", "bitSlipCalibrate",1);
	IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(1000));
	ecc()->writeField("pipeCtrl", "bitSlipCalibrate", 0);
}
//__________________________________________________________________________________________________
/**
 *  Configures an AsAd board's circular buffer.
 */
void EccClient::setAsAdCircularBuffer(const size_t & asadIdx, const uint32_t & startOfMemory, const uint32_t & endOfMemory)
{
	ecc()->connect("device", "ctrl");
	std::string regName = std::string("memStart_asad") + boost::lexical_cast< std::string >(asadIdx);
	ecc()->writeRegister(regName, startOfMemory);
	regName = std::string("memEnd_asad") + boost::lexical_cast< std::string >(asadIdx);
	ecc()->writeRegister(regName, endOfMemory);
	regName = std::string("readPtr_asad") + boost::lexical_cast< std::string >(asadIdx);
	ecc()->writeRegister(regName, startOfMemory);
}
//__________________________________________________________________________________________________
/**
 *  Configures the 4 circular buffers.
 */
void EccClient::setCircularBuffers(const uint32_t & startOfMemory, const uint32_t & endOfMemory)
{
	// Set size of each AsAd board's circular buffer depending on the number of boards
	uint32_t asadBufferSize = endOfMemory - startOfMemory;
	size_t numAsAdBoards = 0;
	for (size_t asadIdx=0; asadIdx <4u; ++asadIdx)
	{
		if (isAsAdBoardActive(asadIdx)) numAsAdBoards++;
	}
	if (3 == numAsAdBoards) numAsAdBoards  = 4;
	asadBufferSize = asadBufferSize/numAsAdBoards;

	ecc()->connect("device", "ctrl");
	try
	{
		uint32_t asadMemStart = startOfMemory;
		for (size_t asadIdx=0; asadIdx <4u; ++asadIdx)
		{
			uint32_t currentBufferSize = (isAsAdBoardActive(asadIdx))?asadBufferSize:0;

			std::string regName = std::string("memStart_asad") + boost::lexical_cast< std::string >(asadIdx);
			ecc()->writeRegister(regName, asadMemStart);
			regName = std::string("memEnd_asad") + boost::lexical_cast< std::string >(asadIdx);
			ecc()->writeRegister(regName, asadMemStart + currentBufferSize);
			regName = std::string("readPtr_asad") + boost::lexical_cast< std::string >(asadIdx);
			ecc()->writeRegister(regName, asadMemStart);

			asadMemStart += currentBufferSize;
		}
	}
	catch (const CmdException &)
	{
		ecc()->writeRegister("memStart", startOfMemory);
		ecc()->writeRegister("memEnd", endOfMemory);
	}
}
//__________________________________________________________________________________________________
/**
 *  Configures parameters common to all 4 circular buffers.
 */
void EccClient::setCircularBufferCommonParameters(const uint32_t & eventsPerInterrupt, const uint32_t & levelAlmostFullMemory,
		const uint32_t & levelAlmostEmptyMemory)
{
	try
	{
		ecc()->writeRegister("interruptRate", eventsPerInterrupt);
	}
	catch (const CmdException &)
	{
	}
	ecc()->writeRegister("aFullLevel", levelAlmostFullMemory);
	ecc()->writeRegister("aEmptyLevel", levelAlmostEmptyMemory);
}
//__________________________________________________________________________________________________
/**
 *  Selects the trigger mode.
 *  It controls the behavior of the CoBo trigger core
 */
void EccClient::setTriggerMode(const std::string & trigger)
{
	/* First, all trigger modes are disabled */
	ecc()->connect("device", "ctrl");
	ecc()->writeField("triggerMask", "external", 0);
	ecc()->writeField("triggerMask", "multiplicity", 0);
	ecc()->writeField("triggerMask", "timing", 0);
	ecc()->writeField("triggerMask", "manual", 0);
	ecc()->writeField("triggerMask", "slowControl", 0);
	ecc()->writeField("triggerMask", "external_2p", 0);
	ecc()->writeField("triggerMask", "multiplicity_2p", 0);
	ecc()->writeField("triggerMask", "timing_2p", 0);
	ecc()->writeField("triggerMask", "manual_2p", 0);
	ecc()->writeField("triggerMask", "slowControl_2p", 0);

	LOG_INFO() << "Selecting '" << trigger << "' primary trigger mode";

	/* the requested trigger mode is enabled */
	if (trigger == "onRequest")
	{
		ecc()->writeField("triggerMask", "manual", 1);
	}
	else if (trigger == "onMultiplicity")
	{
		ecc()->writeField("triggerMask", "multiplicity", 1);
	}
	else if (trigger == "onSlowControl")
	{
		ecc()->writeField("triggerMask", "slowControl", 1);
	}
	else if (trigger == "periodically")
	{
		ecc()->writeField("triggerMask", "timing", 1);
	}
	else if (trigger == "external")
	{
		ecc()->writeField("triggerMask", "external", 1);
	}
	else if (trigger == "onRequest_2p")
	{
		ecc()->writeField("triggerMask", "manual", 1);
		ecc()->writeField("triggerMask", "manual_2p", 1);
	}
	else if (trigger == "onMultiplicity_2p")
	{
		ecc()->writeField("triggerMask", "multiplicity", 1);
		ecc()->writeField("triggerMask", "multiplicity_2p", 1);
	}
	else if (trigger == "onSlowControl_2p")
	{

		ecc()->writeField("triggerMask", "slowControl_2p", 1);
		ecc()->writeField("triggerMask", "slowControl", 1);
	}
	else if (trigger == "periodically_2p")
	{
		ecc()->writeField("triggerMask", "timing", 1);
		ecc()->writeField("triggerMask", "timing_2p", 1);
	}
	else if (trigger == "external_2p")
	{
		ecc()->writeField("triggerMask", "external", 1);
		ecc()->writeField("triggerMask", "external_2p", 1);
	}
}
//__________________________________________________________________________________________________
/**
 *  Sets the mode of the secondary trigger.
 *  It controls the behavior of the CoBo external trigger loggic.
 *  When it triggers, it sends a pulse on connector J6, pin 4 but does not trigger the event readout.
 *  The trigger block uses the same parameters as the primary trigger, but has a different trigger mask.
 */
void EccClient::setSecondaryTriggerMode(const std::string & mode)
{
	/* First, all trigger modes are disabled */
	ecc()->connect("device", "ctrl");
	try
	{
		ecc()->writeField("triggerMask", "external_out", 0);
		ecc()->writeField("triggerMask", "multiplicity_out", 0);
		ecc()->writeField("triggerMask", "timing_out", 0);
		ecc()->writeField("triggerMask", "manual_out", 0);
		ecc()->writeField("triggerMask", "slowControl_out", 0);
	}
	catch (const CmdException & e)
	{
		LOG_WARN() << "Error setting secondary trigger mode: " << e.reason;
		return;
	}

	LOG_INFO() << "Selecting '" << mode << "' secondary trigger mode";

	/* the requested trigger mode is enabled */
	if (mode == "onRequest")
	{
		ecc()->writeField("triggerMask", "manual_out", 1);
	}
	else if (mode == "onMultiplicity")
	{
		ecc()->writeField("triggerMask", "multiplicity_out", 1);
	}
	else if (mode == "onSlowControl")
	{
		ecc()->writeField("triggerMask", "slowControl_out", 1);
	}
	else if (mode == "periodically")
	{
		ecc()->writeField("triggerMask", "timing_out", 1);
	}
	else if (mode == "external")
	{
		ecc()->writeField("triggerMask", "external_out", 1);
	}
}
//_________________________________________________________________________________________________
/*
 * Sets the delay from rising edge of SCW signal to start of multiplicity sliding window.
 * @param delay Number of CKR clock cycles.
 */
void EccClient::setScwMultDelay(const uint64_t & delay)
{
	LOG_INFO() << "Setting 'ScwMultDelay' to " << delay << " CKR clock cycle(s)";
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("ScwMultDelay", delay);
}
//_________________________________________________________________________________________________
/**
 * Sets the delay between SCW low and SCR high.
 * @param delay_10ns Delay, in units of 10 ns.
 */
void EccClient::setScwScrDelay(const uint64_t & delay_10ns)
{
	LOG_INFO() << "Setting 'ScwScrDelay' to " << delay_10ns*10 << " ns";
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("ScwScrDelay", delay_10ns);
}
//_________________________________________________________________________________________________
/**
 * Sets the CoBo in partial/full readout mode.
 * All AGET chips within a CoBo should be configured accordingly (see AGET register 1 bit 14).
 */
void EccClient::setFullReadOutMode(bool isAllChannelRead)
{
	LOG_INFO() << "Selecting " << (isAllChannelRead?"full":"partial") << " readout mode";
	ecc()->connect("device", "ctrl");
	if (isAllChannelRead)
		ecc()->writeField("pipeCtrl", "partialReadOut", 0);
	else
		ecc()->writeField("pipeCtrl", "partialReadOut", 1);
}
//_________________________________________________________________________________________________
void EccClient::setTriggerDelay_10ns(const uint32_t & triggerDelay_10ns)
{
	LOG_INFO() << "Setting CoBo trigger delay to " << triggerDelay_10ns << " x 10 ns";
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("triggerDelay", triggerDelay_10ns);
}
//_________________________________________________________________________________________________
void EccClient::setTriggerDelay2p_10ns(const uint32_t & triggerDelay)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("triggerDelay_2p", triggerDelay);
}
//_________________________________________________________________________________________________
void EccClient::setTriggerTimeOut2p_10ns(int triggerDelay)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("triggerWindowEnd_2p", triggerDelay);
}
//_________________________________________________________________________________________________
void EccClient::setTriggerDeadTime2p_10ns(int triggerDelay)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("triggerWindowStart_2p", triggerDelay);
}
//_________________________________________________________________________________________________
void EccClient::setTriggerTime2p_10ns(int triggerDelay)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("triggerTime_10ns_2p", triggerDelay);
}
//__________________________________________________________________________________________________
/**
 * Configures the subtraction of the moving average and deviation of each AGET's multiplicity signal.
 */
void EccClient::setCoBoMultSubtractMode(bool subtractAverage, const uint8_t & subtractDevs)
{
	ecc()->connect("device", "ctrl");
	try
	{
		LOG_INFO() << "Setting multiplicity subtractAverage to " << subtractAverage << " and subtractDevs to " <<  (short) subtractDevs;
		// When set, the average multiplicity value of the AGET is subtracted from the AGET output to shift the baseline down to 0.
		ecc()->writeField("pipeCtrl", "MultSubtractAverage", subtractAverage?1:0);
		// Twice the number of absolute deviations to subtract from multiplicity signal.  In testing, a value of 4 (2 deviations) worked well to prevent noise from causing triggers.
		ecc()->writeField("pipeCtrl", "MultSubtractDevs", subtractDevs);
	}
	catch (const CmdException & e)
	{
		LOG_WARN() << "Error setting multiplicity subtraction mode: " << e.reason;
		return;
	}
}
//__________________________________________________________________________________________________
/**
 * Configures the 4 bits of the CoBo 'pipeCtrl' register indicating which of the 4 chips to ignore in the multiplicity trigger mode.
 */
void EccClient::setCoBoSuppressMultiplicity(const uint8_t & suppressMultiplicity)
{
	ecc()->connect("device", "ctrl");
	try
	{
		LOG_INFO() << "Setting 'suppressMultiplicity' field to " << std::hex << std::showbase << (unsigned short) suppressMultiplicity << std::noshowbase << std::dec;
		ecc()->writeField("pipeCtrl", "suppressMultiplicity", suppressMultiplicity);
	}
	catch (const CmdException & e)
	{
		LOG_WARN() << "Error setting suppressMultiplicity register field: " << e.reason;
		return;
	}
}
//__________________________________________________________________________________________________
/**
 * Sets 2p mode.
 */
void EccClient::enableMem2pMode(bool valid2pMode)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeField("pipeCtrl", "modeEnable_2p", valid2pMode);
	LOG_INFO() << "Setting 'modeEnable_2p' to " << valid2pMode;
}
//_________________________________________________________________________________________________
/**
 * Configures the value of the 'dataSource' field in the header of the data frames.
 * This Data source ID typically refers to the detector or detector type from which the frame data have been extracted.
 * @param dataSourceId Data source ID.
 */
void EccClient::setDataSourceId(const uint8_t & dataSourceId)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeField("dataSource", "value", dataSourceId);
	LOG_INFO() << "Setting frame data source ID to " << (short) dataSourceId;
}
//_________________________________________________________________________________________________
/**
 * Sets the value of the CoBo index field in the header of the data frames.
 * @param coboId CoBo board identifier.
 */
void EccClient::setCoBoId(const uint8_t & coboId)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeField("COBO_index", "d", coboId);
	LOG_INFO() << "Setting frame CoBo ID to " << (short) coboId;
}
//_________________________________________________________________________________________________
/**
 * Sets the period of the periodic trigger mode.
 * @param period_10ns Period, in units of 10 ns.
 */
void EccClient::triggerPeriod_10ns(const uint32_t & period_10ns)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("triggerPeriod_10ns", period_10ns);
}
//_________________________________________________________________________________________________
/**
 * Launches daqStart (status = true) or daqStop (status = false)
 */
void EccClient::daqCtrl(bool status, bool resetTime)
{
	if (status)
	{
		// Reset timestamp and event counters
		if (resetTime)
		{
			LOG_INFO() << "Resetting timestamp and event counters";
			ecc()->connect("device", "ctrl");
			ecc()->writeRegister("resetTime", 0);
		}

		// Start data sender and interrupt monitor
		LOG_INFO() << "daqStart";
		ecc()->daqStart();

		// Switch to fast control
	    switchToFastControl();
	}
	else
	{
		LOG_INFO() << "daqStop";
		ecc()->stopAsAdPeriodicPulser();
		ecc()->daqStop();
	}
}
//_________________________________________________________________________________________________
void EccClient::daqSetAlwaysFlushData(bool enable)
{
	LOG_INFO() << "Setting data sender 'always flush' policy to " << enable;
	ecc()->setAlwaysFlushData(enable);
}
//_________________________________________________________________________________________________
void EccClient::daqSetCircularBuffersEnabled()
{
	LOG_INFO() << "Selecting circular buffers";
	ecc()->setCircularBuffersEnabled(asadMask_);
}
//_________________________________________________________________________________________________
/**
 * Enables or disables the test mode.
 * In this mode, the frames contain the value of an incremental counter instead of the ADC data.
 */
void EccClient::acqTestMode(bool enable)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeField("pipeCtrl", "testMode", (enable?1:0));
	if (enable)
	{
		MSG_IF(verbose) << "Enabling test acquisition mode!";
	}
}
//__________________________________________________________________________________________________
/**
 * Sets CoBo CKR clock frequency.
 * By design, it must now be set to 25 MHz.
 */
void EccClient::setCoBoReadingClockFrequency(float ckr_MHz)
{
	LOG_INFO() << "Setting CKR clock frequency to " << ckr_MHz << " MHz";
	ecc()->connect("device", "ctrl");
	if (ckr_MHz == 25.0)
	{
		ecc()->writeField("ckrDivider", "d", 2); //  25 MHz
	}
	else if (ckr_MHz == 12.5)
	{
		ecc()->writeField("ckrDivider", "d", 3); //  12.5 MHz
	}
	else if (ckr_MHz == 6.25)
	{
		ecc()->writeField("ckrDivider", "d", 4); //  6.25 MHz
	}
	else
	{
		std::ostringstream msg;
		msg << ckr_MHz << " MHz is not a valid CKR frequency. The CKR frequency must always be set to 25 MHz.";
		throw CmdException(msg.str());
	}
}
//__________________________________________________________________________________________________
/**
 * Sets the manual read data delay.
 * Only has an effect when using a Reduced-CoBo board.
 * Manual delay registers have been disconnected on the micro-TCA CoBo.
 */
void EccClient::setReadDataDelay(int delay)
{
	MSG_IF(verbose) << "Setting manual read data delay to " << delay;
	ecc()->connect("device", "ctrl");
	ecc()->writeField("readDataDelay", "AGET0", delay);
	ecc()->writeField("readDataDelay", "AGET1", delay);
	ecc()->writeField("readDataDelay", "AGET2", delay);
	ecc()->writeField("readDataDelay", "AGET3", delay);
}
//__________________________________________________________________________________________________
/**
 * Enables or disables writing to the hit register.
 */
void EccClient::setWriteHitRegisterEnabled(bool enabled)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeField("pipeCtrl", "writeHitRegisterMaskEnable", enabled?1:0);
}
//_________________________________________________________________________________________________
/**
 * Sets the clocking manager of AsAd for the readout.
 * Applies to all AsAd boards selected for slow control.
 * @see setAsAdSlowControlMask
 */
void EccClient::setAsAdLcm2Register(uint64_t value_LCM2)
{
	MSG_IF(verbose) << "Setting LCM2 register to " << std::hex << std::showbase << value_LCM2 << std::dec << std::noshowbase;
	ecc()->connect("device", "asad");
	ecc()->writeRegister("LCM2", value_LCM2);
	IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(1000));
}
//_________________________________________________________________________________________________
/**
 * Configures the SCA Write clock frequency.
 *
 * @param pllConfigEnabled Whether to configure CoBo PLL device (for debugging purposes).
 */
void EccClient::setWritingClockFrequency(float ckw_MHz, bool pllConfigEnabled, const get::cobo::PllRegisterMap & regs)
{
	LOG_INFO() << "Setting CKW clock frequency to " << ckw_MHz << " MHz";
	if (not pllConfigEnabled)
	{
		LOG_WARN() << "CoBo PLL device will not be configured!";
	}
	ecc()->coboSetWritingClockFrequency(ckw_MHz, pllConfigEnabled, regs);
}
//_________________________________________________________________________________________________
/**
* Sets the clocking manager of AsAd for the sampling.
* Applies to all AsAd boards selected for slow control.
* @see setAsAdSlowControlMask
*/
void EccClient::setAsAdLcm1Register(const uint64_t & value)
{
	MSG_IF(verbose) << "Setting LCM1 register to "  << std::hex << std::showbase << value << std::dec << std::noshowbase;
	ecc()->connect("device", "asad");
	ecc()->writeRegister("LCM1", value);
}
//_________________________________________________________________________________________________
/**
 * Sets function of LEMO connector on CoBo front panel.
 */
void EccClient::setCoBoLemoMode(const size_t lemoIndex, const uint8_t & mode)
{
	//LOG_DEBUG() << "Setting LEMO no. " << lemoIndex << " to mode " << (uint16_t) mode;
	ecc()->connect("device", "ctrl");
	const std::string fieldName = std::string("L") + boost::lexical_cast< std::string >(lemoIndex);
	ecc()->writeField("lemoMux", fieldName, mode);
}
//_________________________________________________________________________________________________
/**
 * Sets the ADC to send a constant vector
 */
void EccClient::adcTest(ModeADC modeADC)
{
	ecc()->connect("device", "asad");
	ecc()->writeRegister("ADC2", 0x2A0);
	ecc()->writeRegister("ADC3", 0x001);
	ecc()->writeRegister("ADC4", 0x000);
	ecc()->writeRegister("ADC5", (modeADC == SDR)?0x403:0x401);
}
//_________________________________________________________________________________________________
/**
 * Checks whether AsAd board is plugged in into the CoBo board.
 * Reads CoBo read-only register giving status of the AsAd boards.
 * @param asadIdx Index of the AsAd board.
 */
bool EccClient::isAsAdConnected(const size_t & asadIdx)
{
	ecc()->connect("device", "ctrl");
	const std::string fieldName = std::string("PLG") + boost::lexical_cast< std::string >(asadIdx);
	bool isConnected = (1 == ecc()->readField("asadConnection", fieldName));
	return isConnected;
}
//_________________________________________________________________________________________________
/**
 * Checks whether AsAd board is connected to its power supply.
 * Reads AsAd alarm bits (value in CoBo reversed w.r.t. AsAd?).
 * @param asadIdx Index of the AsAd board.
 */
bool EccClient::isAsAdPowerSupplied(const size_t & asadIdx)
{
	return getAsAdAlarmFlag(asadIdx);
}
//_________________________________________________________________________________________________
/**
 * Checks whether AsAd board was powered on.
 * @see GET-AS-002-0004
 */
bool EccClient::isAsAdPoweredOn(const size_t & asadIdx)
{
	ecc()->connect("device", "ctrl");
	const std::string fieldName = std::string("powerON") + boost::lexical_cast< std::string >(asadIdx);
	bool pwu = (1 == ecc()->readField("asadEnable", fieldName));
	return pwu;
}
//_________________________________________________________________________________________________
/**
 * Writes a value in a register of the "aget" device.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param registerName Name of the register.
 * @param value Value to write.
 */
void EccClient::agetRegWSC(const std::size_t & asadIdx, const std::size_t & agetIdx, const std::string & registerName, uint64_t value)
{
	try
	{
		// Select AGET to write to
		ecc()->connect("device", "ctrl");
		uint16_t mask = 1 << (4*asadIdx + agetIdx);
		ecc()->writeField("aget_sc_configuration", "SC_En", mask);
		// Write register value
		ecc()->connect("device", "aget");
		ecc()->writeRegister(registerName, value);
	}
	catch (const CmdException & e)
	{
		LOG_ERROR() << "Could not set AGET no. " << asadIdx*4+agetIdx << " register '" << registerName << "': " << e.reason;
	}
	catch (const std::exception& e)
	{
		LOG_ERROR() << "Could not set AGET no. " << asadIdx*4+agetIdx << " register '" << registerName << "': " << e.what();
	}
}
//_________________________________________________________________________________________________
/**
 * Writes a value in a register of all the enabled "AGET" devices.
 * @param registerName Name of the register.
 * @param value Value to write.
 */
void EccClient::agetRegWSC(const std::string & registerName, uint64_t value)
{
	try
	{
		// Write register value
		ecc()->connect("device", "aget");
		ecc()->writeRegister(registerName, value);
	}
	catch (const CmdException & e)
	{
		LOG_ERROR() << "Could not set AGET register '" << registerName << "': " << e.reason;
	}
	catch (const std::exception& e)
	{
		LOG_ERROR() << "Could not set AGET no. register '" << registerName << "': " << e.what();
	}
}
//_________________________________________________________________________________________________
/**
 * Reads a value from a register of the "aget" device.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param registerName Name of the register.
 * @return Returns the register value read.
 */
uint64_t EccClient::agetRegRSC(const std::size_t & asadIdx, const std::size_t & agetIdx, const std::string & registerName)
{
	uint64_t regValue = UINT64_C(0);
	try
	{
		// Select AGET to read from
		ecc()->connect("device", "ctrl");
		// FIXME: also set aget_sc_configuration/SC_En mask ?
		ecc()->writeField("agetChipSelect", "d", 4*asadIdx + agetIdx);
		// Read register value
		ecc()->connect("device", "aget");
		regValue = ecc()->readRegister(registerName);
	}
	catch (const CmdException & e)
	{
		LOG_ERROR() << "Could not read AGET no. " << asadIdx*4+agetIdx << " register '" << registerName << "': " << e.reason;
	}
	catch (const std::exception& e)
	{
		LOG_ERROR() << "Could not read AGET no. " << asadIdx*4+agetIdx << " register '" << registerName << "': " << e.what();
	}
	return regValue;
}
//_________________________________________________________________________________________________
/**
 * Sets ReadIfHit registers for a given AGET chip.
 * @param asadIdx Index of the AsAd board the AGET belongs to (0 to 3).
 * @param agetIdx Index of the AGET chip within its AsAd board (0 to 3).
 * @param val1_4 Value of the readIfHitMask1_4 register.
 * @param val5_36 Value of the readIfHitMask5_36 register.
 * @param val37_68 Value of the readIfHitMask37_68 register.
 */
void EccClient::setReadIfHitMask(const std::size_t & asadIdx, const std::size_t & agetIdx, const uint32_t & val1_4, const uint32_t &  val5_36, const uint32_t &  val37_68)
{
	ecc()->connect("device", "ctrl");
	// Select AGET chip
	ecc()->writeRegister("maskSelect", long(asadIdx*4+agetIdx));
	// Write register values
	LOG_INFO() << "Setting 'readIfHitMask' registers for AGET " << asadIdx*4+agetIdx << " to " << std::showbase << std::hex << val1_4 << ", " << val5_36 << ", " << val37_68 << std::dec;
	ecc()->writeRegister("readIfHitMask1_4", val1_4);
	ecc()->writeRegister("readIfHitMask5_36", val5_36);
	ecc()->writeRegister("readIfHitMask37_68", val37_68);
}
//__________________________________________________________________________________________________
/**
 * Sets ReadAlways registers for a given AGET chip.
 * @param asadIdx Index of the AsAd board the AGET belongs to (0 to 3).
 * @param agetIdx Index of the AGET chip within its AsAd board (0 to 3).
 * @param val1_4 Value of the readAlwaysMask1_4 register.
 * @param val5_36 Value of the readAlwaysMask5_36 register.
 * @param val37_68 Value of the readAlwaysMask37_68 register.
 */
void EccClient::setReadAlwaysMask(const std::size_t & asadIdx, const std::size_t & agetIdx, const uint32_t & val1_4, const uint32_t & val5_36, const uint32_t & val37_68)
{
	ecc()->connect("device", "ctrl");
	// Select AGET chip
	ecc()->writeRegister("maskSelect", asadIdx*4+agetIdx);
	// Write register values
	LOG_INFO() << "Setting 'readAlwaysMask' registers for AGET " << asadIdx*4+agetIdx << " to " << std::showbase << std::hex << val1_4 << ", " << val5_36 << ", " << val37_68 << std::dec;
	ecc()->writeRegister("readAlwaysMask1_4", val1_4);
	ecc()->writeRegister("readAlwaysMask5_36", val5_36);
	ecc()->writeRegister("readAlwaysMask37_68", val37_68);
}
//__________________________________________________________________________________________________
/**
 * blinking of one led
 */
void EccClient::led(void)
{
	ecc()->connect("device", "debug");
	ecc()->writeField("led", "l0", 1);
	sleep_ms(400);
	ecc()->writeField("led", "l0", 0);
	sleep_ms(400);
}
//__________________________________________________________________________________________________
/**
 * Switches to fast control.
 * This function switches the state machine of the CoBo AGET controller
 * and the state machine on the AGET to the fast control mode,
 * i.e. the reading/writing of the hit register
 */
void EccClient::switchToFastControl()
{
	// Resets the state machine for the AGET slow control
	ecc()->connect("device", "ctrl");
	ecc()->writeField("agetController", "reset", 1); // Will reset the error flags and the busy flag (needs to be pulsed: 0,1,0).
	ecc()->writeField("agetController", "reset", 0);

	// Makes a fake slow control to send the three pulses
	// necessary for AGET to switch to the fast control mode
	ecc()->writeField("aget_sc_configuration", "SC_En", 0); // Mask indicating which of the 16 AGet's the slow control commands should be broadcast to.
	ecc()->writeField("aget_sc_configuration", "pulseEnable", agetMask_); // Needed when switching between slow control mode and hit register mode. Bits written into the first Byte read by the FPGA.
	ecc()->writeField("agetController", "hitRegister", 1); // If set, switches to hit register mode.
	ecc()->writeField("agetController", "calibration", 0); // If set, calibration delays are automatically computed.

	ecc()->writeRegister("aget_sc_rw", 1); // Triggers the read (1) or write (0) operation of the AGET register.

	// enables the control of the AGET and disables the pulse generation
	// No 3-pulse should be sent during a fast control of the AGET
	ecc()->writeField("aget_sc_configuration", "SC_En", agetMask_);
	ecc()->writeField("aget_sc_configuration", "pulseEnable", 0);
}
//__________________________________________________________________________________________________
/**
 * Switches to slow control.
 * This function switches the state machine of the CoBo AGET controller
 * and the state machine of the AGET to the slow control mode
 */
void EccClient::switchToSlowControl()
{
	/* resets the state machine for the AGET slow control
	 * and disables the firmware hand on the AGET slow control */
	ecc()->connect("device", "ctrl");
	ecc()->writeField("pipeCtrl", "firmwareEnable", 0);
	//daqCtrl(false);
	ecc()->writeField("agetController", "reset", 1);
	ecc()->writeField("agetController", "reset", 0);

	/* makes a fake slow control to send the three pulses
	 * necessary for AGET to switch to the slow control mode
	 */
	ecc()->writeField("aget_sc_configuration", "SC_En", 0);
	ecc()->writeField("aget_sc_configuration", "pulseEnable", agetMask_);
	ecc()->writeField("agetController", "hitRegister", 0);
	ecc()->writeField("agetController", "calibration", 0);

	ecc()->writeRegister("aget_sc_rw", 1);

	/* enables the control of the AGET and disables the pulse generation
	 * No 3-pulse should be sent during a fast control of the AGET
	 */
	ecc()->writeField("aget_sc_configuration", "SC_En", agetMask_);
	ecc()->writeField("aget_sc_configuration", "pulseEnable", 0);
}
//__________________________________________________________________________________________________
/**
 * Injects charge in mV.
 * @param deltaV_mV the charge to inject is expressed as a delta voltage in mV
 */
void EccClient::injectCharge(long deltaV_mV)
{
	try
	{
		LOG_INFO() << "injectCharge(" << std::dec << deltaV_mV << " mV)";
		ecc()->triggerAsAdPulser(deltaV_mV);
	}
	catch (const CmdException & e)
	{
		LOG_ERROR() << e.what() << ": " << e.reason;
	}
	catch (const Ice::ConnectionRefusedException&)
	{
		LOG_ERROR() << "Connection refused -> Is ECC server down? Try restarting ECC server";
	}
	catch (const std::exception& e)
	{
		LOG_ERROR() << e.what();

	}
	catch (const std::string& msg)
	{
		LOG_ERROR() << msg;
	}
	catch (const char* msg)
	{
		LOG_ERROR() << msg;
	}
	catch (...)
	{
		LOG_ERROR() << "Unknown error";
	}
}
//_________________________________________________________________________________________________
/**
 * Enables the systematic readout of the FPN channels in both CoBo and the currently selected AGET chips.
 */
void EccClient::enableFPNRead(bool isFPNRead)
{
	MSG_IF(verbose) << (isFPNRead?"Enabling":"Disabling") << " readout of FPN channels";
	ecc()->connect("device", "ctrl");
	ecc()->writeField("pipeCtrl", "isFPNRead", isFPNRead?1:0);
	ecc()->connect("device", "aget");
	ecc()->writeField("reg1", "FPNreadout", isFPNRead?1:0);
}
//_________________________________________________________________________________________________
/**
 * Triggers immediately the reading of the AGET.
 * It is done by writing whatever value of the register "manualTriggerNow"
 * The trigger has to be in "on Request" mode
 */
void EccClient::triggerNow()
{
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("manualTriggerNow", 1);
}
//_________________________________________________________________________________________________
/**
 * Sleeps for time_ms, in milliseconds
 * It uses the Ice sleep function
 */
void EccClient::sleep_ms(int time_ms)
{
	IceUtil::ThreadControl::sleep(IceUtil::Time::milliSeconds(time_ms));
}
//_________________________________________________________________________________________________
/**
 * Set the multiplicity threshold in CoBo register
 */
void EccClient::setMultThreshold(int thres)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("multThreshold", thres);
}
//_________________________________________________________________________________________________
/**
 * Set the multiplicity threshold in CoBo register
 */
void EccClient::setMultThreshold2p(int thres){
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("multThreshold_2p", thres);
}
//_________________________________________________________________________________________________
/**
 * Sets width of multiplicity sliding window.
 */
void EccClient::setSlidingWin(int slid){
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("multWindSize", slid);
}
//__________________________________________________________________________________________________
/**
 * Sets the CoBo 'readDepth' register for the number of time buckets to readout (10 bits)
 */
void EccClient::setCoBoReadDepth(const uint16_t & depth)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeRegister("readDepth", depth);
}
//__________________________________________________________________________________________________
/**
 * Sets the CoBo 'readOffset' register for the time bucket to readout from (9 bits)
 */
void EccClient::setCoBoReadOffset(const uint16_t & offset)
{
	ecc()->connect("device", "ctrl");
	try
	{
		ecc()->writeField("readOffset", "value", offset);
	}
	catch (const CmdException & e)
	{
		LOG_WARN() << e.what(); // This register was only described on Nov 25th, 2013
	}
}
//__________________________________________________________________________________________________
/**
 * Sets the CoBo 'zeroSuppressionEnabled' register to enable zero suppression.
 * @param enabled Whether to enable or disable zero suppression.
 * @param inverted  Whether to suppress values below (false) or above thresholds.
 */
void EccClient::setCoBoZeroSuppressionMode(const bool & enabled, const bool & inverted)
{
	ecc()->connect("device", "ctrl");
	ecc()->writeField("pipeCtrl", "zeroSuppressionEnable", enabled?1:0);
	try
	{
		ecc()->writeField("pipeCtrl", "zeroSuppressInvert", inverted?1:0);
	}
	catch (const CmdException & e)
	{
		LOG_WARN() << e.what(); // This bit was described after email from N. U. on May 20th, 2014.
	}
}
//__________________________________________________________________________________________________
/**
 * Sets the threshold for zero suppression for the given channel.
 * @param asadIdx AsAd ID
 * @param agetIdx AGET chip ID
 * @param channelIdx Channel ID
 * @param threshold Threshold below which to suppress data.
 */
void EccClient::setCoBoZeroSuppressionThreshold(const size_t & asadIdx, const size_t & agetIdx, const size_t channelIdx, const uint16_t & threshold)
{
	std::ostringstream deviceName;
	deviceName << "zeroSuppress";
	if (asadIdx > 0)
	{
		deviceName << "_AsAd" << asadIdx;
	}
	ecc()->connect("device", deviceName.str());
	std::ostringstream regName;
	regName << "aget" << agetIdx << "ch" << channelIdx;
	ecc()->writeField(regName.str(), "threshold", threshold);
}
//__________________________________________________________________________________________________
/**
 * Configures the AsAd TCM device.
 *   - selects step voltage or current pulse mode
 *   - sets capacitance
 *   - selects trigger mode (external or slow control)
 *   - sets triggering delay
 *
 * @param TC  Selects test mode (TC=0, step voltage) or calibration mode (TC=1, current pulse).
 * @param RNG Selects 1 pF (RNG=0) or 11 pF external capacitor (for TC=1).
 * @param TRG Selects slow control (0) or external (1) triggering mode of AsAd voltage generator.
 * @param defaultVoltage_mV Default output voltage (w.r.t. Vref), in mV.
 * @param voltage_mV Output voltage (w.r.t. Vref), in mV.
 *
 * @see Section 1.4 of GET-AS-002-0005: "AsAd Test and Calibration Management".
 */
void EccClient::asadConfigureTcmDevice(const bool & TC, const bool & RNG, const bool & TRG, const uint8_t & DLY,
		const float & defaultVoltage_mV, const float & voltage_mV)
{
	MSG_IF(verbose) << "Configuring AsAd generator: mode T/C=" << TC << " RNG=" << RNG << " TRG=" << TRG << " DLY=" << (uint16_t) DLY;

	ecc()->connect("device", "asad");
	ecc()->writeField("TCM0", "Calib", TC);
	ecc()->writeField("TCM0", "Rg", RNG);
	ecc()->writeField("TCM0", "Delay", DLY);

	// Enable AsAd TCM device slow control triggering mode
	ecc()->writeField("TCM0", "Trg", 0);

	// Set default output voltage
	MSG_IF(verbose) << "Setting AsAd generator default output voltage to " <<  (1102.0f + defaultVoltage_mV) << " mV";
	ecc()->setAsAdPulserDefaultVoltage(defaultVoltage_mV);

	// Enable AsAd TCM device external triggering mode if needed
	if (TRG)
	{
		MSG_IF(verbose) << "Setting AsAd generator output voltage to " <<  (1102.0f + voltage_mV) << " mV";
		ecc()->connect("device", "asad");
		ecc()->writeField("TCM0", "Trg", TRG);
		ecc()->triggerAsAdPulser(voltage_mV);
	}
	try
	{
		ecc()->configureAsAdExternalPulser(TRG, voltage_mV);
	}
	catch (const ::Ice::OperationNotExistException & e) // This operation was added later on
	{
		LOG_WARN() << e.what();
	}
}
//__________________________________________________________________________________________________
void EccClient::asadMonitoringReset()
{
	ecc()->connect("device", "asad");

	const uint16_t status = getAsAdAlertStatus();
	DBG_IF(verbose) << "AsAd monitoring alert status: " << std::hex << std::showbase << status
			<< " (" << get::asad::utl::stringifyMonitoringAlertStatus(status) << ')' << std::dec;

	// Reset all registers to default value
	ecc()->writeRegister("configuration2", 0x80);
}
//__________________________________________________________________________________________________
/**
 * Disable all criteria for entering auto-protection.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring.
 */
void EccClient::asadDisableAutoProtection(const size_t & asadIdx)
{
	const uint8_t interruptMask1 = 0xFF;
	LOG_DEBUG() << "Setting AsAd no. " << asadIdx << " auto-protection level to " << std::hex << std::showbase << (uint16_t) interruptMask1 << std::dec;
	ecc()->connect("device", "asad");
	ecc()->writeRegister("interruptMask1", interruptMask1);
}
//__________________________________________________________________________________________________
/**
 * Defines AsAd auto-protection level.
 * See GET-AS-002-0004, section 2.1.2
 */
void EccClient::asadSetAutoProtectionLevel(const uint8_t & interruptMask1, const uint8_t & /* interruptMask2 */)
{
	ecc()->connect("device", "asad");
	// Set auto-protection level
	DBG_IF(verbose) << "Setting AsAd auto-protection level to " << std::hex << std::showbase << (uint16_t) interruptMask1 << std::dec;
	ecc()->writeRegister("interruptMask1", interruptMask1);
	//ecc()->writeRegister("interruptMask2", interruptMask2); // Should be left at default value
}
//__________________________________________________________________________________________________
/**
 * Defines operational mode for AsAd monitoring.
 * See GET-AS-002-0004, section 2.1.3
 *   - Defines auto-protection level
 *   - Defines operational mode for monitoring
 */
void EccClient::asadSetMonitoringMode(const uint8_t & config1, const uint8_t & config2, const uint8_t & config3)
{
	ecc()->connect("device", "asad");

	// ecc()->writeRegister("intTempOffset", 0x0);
	// ecc()->writeRegister("extTempOffset", 0x0);

	// Define operational mode
	DBG_IF(verbose) << "Setting AsAd register 'configuration3' to " << std::hex << std::showbase << (uint16_t) config3 << std::dec;
	ecc()->writeRegister("configuration3", config3);

	// configuration2:
	//  0x0: Round-robin with averaging
	// 0x34: Single channel without averaging
	// 0x20: Round-robin without averaging
	DBG_IF(verbose) << "Setting AsAd register 'configuration2' to " << std::hex << std::showbase << (uint16_t) config2 << std::dec;
	ecc()->writeRegister("configuration2", config2);

	DBG_IF(verbose) << "Setting AsAd register 'configuration1' to " << std::hex << std::showbase << (uint16_t) config1 << std::dec;
	ecc()->writeRegister("configuration1", config1);
}
//__________________________________________________________________________________________________
void EccClient::startMonitorAlarm()
{
	ecc()->connect("device", "asad");
	ecc()->writeRegister("interruptMask1", 0xA3); //0x23); //0x30);
	sleep_ms(10);
}
//__________________________________________________________________________________________________
/**
 * Checks AsAd AL bit.
 * @param asadIdx AsAd board Identifier.
 */
bool EccClient::getAsAdAlarmFlag(const size_t & asadIdx)
{
	ecc()->connect("device", "ctrl");
	const uint8_t alarmBits = ecc()->readField("asadStatus", "alarm");
	const bool alarmBitSet = utl::BitFieldHelper< uint8_t >::getBit(alarmBits, asadIdx);
	return alarmBitSet;
}
//__________________________________________________________________________________________________
/**
 * Reads AsAd monitoring interrupt status register "interruptStatus1".
 * Needed for reset after an alert has been triggered.
 * @return Value of the interruptStatus1 register before reset.
 */
uint8_t EccClient::getAsAdAlertStatus()
{
	ecc()->connect("device", "asad");
	return (uint8_t) ecc()->readRegister("interruptStatus1");
}
//__________________________________________________________________________________________________
/**
 * Stops AsAd monitoring.
 */
void EccClient::asadMonitoringStop()
{
	ecc()->connect("device", "asad");
	ecc()->writeRegister("configuration1", 0x2c);
}
//__________________________________________________________________________________________________
/**
 * Reads AsAd VDD monitoring register.
 * @return VDD voltage, in straight binary format.
 * See AsAd Power Supply and Monitoring document.
 */
uint8_t EccClient::readVdd(const std::size_t & asadIdx)
{
	try
	{
		// Select AsAd board
		ecc()->connect("device", "ctrl");
		ecc()->writeField("asadSelect", "d", asadIdx);
		// Read register
		ecc()->connect("device", "asad");
		uint8_t val = (uint8_t) ecc()->readRegister("VDD");
		DBG_IF(verbose) << "AsAd " << asadIdx << " VDD voltage: "  << std::showbase << std::hex << (uint16_t) val << std::dec;
		return val;
	}
	catch (...)
	{
		LOG_ERROR() << "Error reading VDD (AsAd #" << asadIdx << ") register";
	}
	return 0;
}
//__________________________________________________________________________________________________
/**
 * Reads IDD AsAd monitoring register.
 */
uint8_t EccClient::readIdd(const std::size_t & asadIdx)
{
	try
	{
		// Select AsAd board
		ecc()->connect("device", "ctrl");
		ecc()->writeField("asadSelect", "d", asadIdx);
		// Read register
		ecc()->connect("device", "asad");
		uint8_t val = (uint8_t) ecc()->readRegister("IDD");
		DBG_IF(verbose) << "AsAd " << asadIdx << " IDD current: "  << std::showbase << std::hex << (uint16_t) val << std::dec;
		return val;
	}
	catch (CmdException & e)
	{
		LOG_ERROR() << e.reason;
	}
	catch (...)
	{
		LOG_ERROR() << "Error reading IDD (Asad) register";
	}
	return 0;
}
//__________________________________________________________________________________________________
/**
 * Reads VDD1 AsAd monitoring register.
 */
uint8_t EccClient::readVdd1(const std::size_t & asadIdx)
{
	try
	{
		// Select AsAd board
		ecc()->connect("device", "ctrl");
		ecc()->writeField("asadSelect", "d", asadIdx);
		// Read register
		ecc()->connect("device", "asad");
		uint8_t val = (uint8_t) ecc()->readRegister("VDD1");
		DBG_IF(verbose) << "AsAd " << asadIdx << " VDD1 voltage: "  << std::showbase << std::hex << (uint16_t) val << std::dec;
		return val;
	}
	catch (CmdException & e)
	{
		LOG_ERROR() << e.reason;
	}
	catch (...)
	{
		LOG_ERROR() << "Error reading VDD1 (AsAd #" << asadIdx << ") register";
	}
	return 0;
}
//__________________________________________________________________________________________________
/**
 * Reads internal temperature AsAd monitoring register.
 */
int8_t EccClient::readIntTemp(const std::size_t & asadIdx)
{
	try
	{
		// Select AsAd board
		ecc()->connect("device", "ctrl");
		ecc()->writeField("asadSelect", "d", asadIdx);
		// Read register
		ecc()->connect("device", "asad");
		int8_t val = (int8_t) ecc()->readRegister("internalTemperature");
		//DBG_IF(verbose) << "AsAd " << asadIdx << " internal temperature: " << std::showbase << std::hex << (int16_t) val << std::dec;
		return val;
	}
	catch (CmdException & e)
	{
		LOG_ERROR() << e.reason;
	}
	catch (...)
	{
		LOG_ERROR() << "Error reading Internal Temperature (AsAd #" << asadIdx << ") register";
	}
	return 0;
}
//__________________________________________________________________________________________________
/**
 * Reads external temperature AsAd monitoring register.
 */
int8_t EccClient::readExtTemp(const std::size_t & asadIdx)
{
	try
	{
		// Select AsAd board
		ecc()->connect("device", "ctrl");
		ecc()->writeField("asadSelect", "d", asadIdx);
		// Read register
		ecc()->connect("device", "asad");
		int8_t val = (int8_t) ecc()->readRegister("externalTemperature");
		//DBG_IF(verbose) << "AsAd " << asadIdx << " external temperature: "  << std::showbase << std::hex << (int16_t) val << std::dec;
		return val;
	}
	catch (CmdException & e)
	{
		LOG_ERROR() << e.reason;
	}
	catch (...)
	{
		LOG_ERROR() << "Error reading External Temperature (AsAd #" << asadIdx << ") register";
	}
	return 0;
}
//__________________________________________________________________________________________________
/**
 * Sets the ISPM multiplexer to debug the signals on AsAd.
 * See Fig. 14 and Tables 14-20 page 17 of "AsAd Firmware Description".
 * See "AsAd External Instruments Interface".
 * @param moduleName Name of device.
 * @param signalNameISP1 Name of ISP1 signal.
 * @param signalNameISP2 Name of ISP2 signal.
*/
void EccClient::asadSetISPM(const std::size_t & /* coboIdx */, const std::size_t & /* asadIdx */, const std::string & moduleName, const std::string & isp1SignalName, const std::string & isp2SignalName)
{
	// Lookup register value.
	const uint64_t value = get::asad::utl::buildIspmRegisterValue(moduleName, isp1SignalName, isp2SignalName);

	// Set register
	MSG_IF(verbose) << "Setting ISPM register to " << std::hex << std::showbase << value << std::dec << std::noshowbase;
	ecc()->connect("device", "asad");
	ecc()->writeRegister("ISPM", value);
}
//__________________________________________________________________________________________________
/**
 *   Monitoring limitation
 */
void EccClient::wVddLimit(uint8_t low, uint8_t high)
{
	try
	{
		ecc()->connect("device", "asad");
		ecc()->writeRegister("VDDLowLimit", low);
		ecc()->writeRegister("VDDHighLimit", high);
	}
	catch (CmdException & e)
	{
		LOG_ERROR() << "Error loading VDD limits (AsAd) registers: " << e.reason;
	}
}
//__________________________________________________________________________________________________
void EccClient::wIddLimit(uint8_t low, uint8_t high)
{
	try
	{
		ecc()->connect("device", "asad");
		ecc()->writeRegister("IDDLowLimit", low);
		ecc()->writeRegister("IDDHighLimit", high);
	}
	catch (CmdException & e)
	{
		LOG_ERROR() << "Error loading IDD limits (AsAd) registers: ";
		LOG_ERROR() << e.reason;
	}
}
//__________________________________________________________________________________________________
void EccClient::wVdd1Limit(uint8_t low, uint8_t high)
{
	try
	{
		ecc()->connect("device", "asad");
		ecc()->writeRegister("VDD1LowLimit", low);
		ecc()->writeRegister("VDD1HighLimit", high);
	}
	catch (CmdException & e)
	{
		LOG_ERROR() << "Error loading VDD1 limits (AsAd) registers: ";
		LOG_ERROR() << e.reason;
	}
}
//__________________________________________________________________________________________________
void EccClient::wIntTempLimit(uint8_t low, uint8_t high)
{
	try
	{
		ecc()->connect("device", "asad");
		ecc()->writeRegister("intTempLowLimit", low);
		ecc()->writeRegister("intTempHighLimit", high);
	}
	catch (CmdException & e)
	{
		LOG_ERROR() << "Error loading internal temperature limits (AsAd) registers: ";
		LOG_ERROR() << e.reason;
	}
}
//__________________________________________________________________________________________________
void EccClient::wExtTempLimit(uint8_t low, uint8_t high)
{
	try
	{
		ecc()->connect("device", "asad");
		ecc()->writeRegister("extTempLowLimit", low);
		ecc()->writeRegister("extTempHighLimit", high);
	}
	catch (CmdException & e)
	{
		LOG_ERROR() << "Error loading external temperature limits (AsAd) registers: ";
		LOG_ERROR() << e.reason;
	}
}
//__________________________________________________________________________________________________
