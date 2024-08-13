/* =====================================================================================================================
 * TestController.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class TestController
 * Created on: 9 sept. 2010 at Irfu/Sedi/Ldef, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Frederic Druillole (frederic.druillole@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GET (General Electronics for TPC) data acquisition framework @ CEA-Irfu.
 * Find out more on: http://www-actar-get.cea.fr
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
#include "get/rc/CoBoConfigParser.h"
#include "get/TestController.h"
#include "get/DefaultPortNums.h"
#include "get/GetControllerException.h"
#include "get/EccClient.h"
#include "get/CatchAllExceptions.h"
#include "get/ControllerOptions.h"
#include "get/WorkspaceManager.h"
#include "get/daq/FrameStorage.h"
#include "get/utl/Time.h"

#include "get/asad/utl.h"
#include "mdaq/utl/net/SocketAddress.h"
#include "CCfg/Attribute.h"
#include "mfm/Serializer.h"
#include "mfm/Frame.h"
#include "mfm/Field.h"
#include "mfm/BitField.h"
#include "mfm/Exception.h"
#include "mfm/FrameFormat.h"
#include "mdaq/utl/Server.h"
#include "utl/BitFieldHelper.hpp"
#include "utl/Logging.h"

#include <boost/cstdint.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/case_conv.hpp>
namespace ba = boost::algorithm;
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
namespace io = boost::iostreams;

#include <QDateTime>
#include <QTextStream>
#include <QTimer>
#include <QMutexLocker>
#include <QDir>
#include <QUuid>

using namespace get;
using namespace ::mdaq::utl::net;

#include <iostream>
#include <sstream>
#include <ctime>
using std::size_t;

#define DBG_IF(cond)  if(!(cond)) ; else LOG_DEBUG()
#undef MSG_IF
#define MSG_IF(cond)  DBG_IF(cond)

/**
 * Static variables
 */
bool TestController::verbose = true;

inline uint64_t logtwo(uint64_t seed)
{
	uint64_t result = 0;
	int maxBits = 64;
	for (int i=0; i < maxBits; ++i)
	{
		switch (seed >> result)
		{
		case 1:
			return result;

		case 0:
			result -= maxBits >> i;
			break;

		default:
			result += maxBits >> i;
			break;
		}
	}
	return result;
}

/**
 * Constructs a controller for a test.
 * @param testId Identifier of the test.
 */
TestController::TestController(const std::string & testName, DataBloc& dataBloc, QObject* parentObject, IFrameProcessor* frameDisplay) :
	QThread(parentObject),
	stopLoop(false),
	nbEvent(0),
	dataBloc(dataBloc),
	coboInstanceName_("CoBo"),
	isRunFinished_(true),
	isSaveRawFrameEnabled_(true),
	isDataMonitoringEnabled_(true),
	frameDisplay(frameDisplay)
{
	LOG_DEBUG() << "Creating test controller";

	// Store test parameters
	testSummaryPtr = WorkspaceManager::instance().readTest(testName);

	// Load run conditions
	CCfg::CConfig cfg(&runCfgDoc_.loadFromFile(testSummary().paramConfigPath()));
	if (cfg("Node", "CoBo").exists())
	{
		runCfg_.setBaseAttribute(&cfg("Node", "CoBo").find());
		coboInstanceName_ = "Instance";
	}
	else
	{
		runCfg_.setBaseAttribute(cfg.getBaseAttribute());
		coboInstanceName_ = "CoBo";
	}
	try
	{
		// Setting verbose mode from parameter file
		setVerbose(runCfg_("ECC")("Verbose").getBoolValue());
	}
	catch (const CCfg::CConfig::NotFound&)
	{ /* keep default value*/
	}

	// See which AsAd/AGET are active
	setMasks(coboConfig(0));

	MSG_IF(verbose) << "ECC -> " << testSummary().eccEndpoint();
	MSG_IF(verbose) << "TARGET -> " << testSummary().targetEndpoint();

	// Enable data display or not
	bool isDisplayEnabled = true;
	try
	{
		isDisplayEnabled = runCfg_("ECC")("DataMonitor").getBoolValue();
	}
	catch (const CCfg::CConfig::NotFound &) {}
	ControllerOptions::instance().setDataDisplayEnabled(isDisplayEnabled);
	MSG_IF(isDisplayEnabled) << "Data display is enabled";

	LOG_DEBUG() << "Verbose mode is " << (isVerbose()?"":"NOT ") << "enabled";
	LOG_DEBUG() << "Data display is " << (ControllerOptions::instance().isDataDisplayEnabled()?"":"NOT ") << "enabled";
}
//__________________________________________________________________________________________________
/**
 * Destructor.
 */
TestController::~TestController()
{
	LOG_DEBUG() << "Destroying test controller";
}
//__________________________________________________________________________________________________
/**
 * Returns mask of enabled AsAd board(s).
 */
uint16_t TestController::asadMask() const
{
	return asadMask_;
}
//__________________________________________________________________________________________________
/**
 * Checks in configuration which AsAd boards and AGET chips to enable.
 * @param coboCfg Configuration of the CoBo board.
 */
void TestController::setMasks(CCfg::CConfig coboCfg)
{
	asadMask_ = 0;
	agetMask_ = 0;

	if (coboCfg("isActive").getBoolValue())
	{
		// Loop over 4 AsAd boards
		for (size_t asadIdx=0; asadIdx < 4; ++asadIdx)
		{
			try
			{
				bool isBoardActive = coboCfg("AsAd", asadIdx)("Control")("isActive").getBoolValue();
				if (isBoardActive)
				{
					::utl::BitFieldHelper<uint16_t>::setField(asadMask_, asadIdx, 1, 1);
					// Loop over 4 AGET chips
					for (size_t agetIdx=0; agetIdx < 4; ++agetIdx)
					{
						bool isChipActive = coboCfg("AsAd", asadIdx)("Aget", agetIdx)("Control")("isActive").getBoolValue();
						::utl::BitFieldHelper<uint16_t>::setField(agetMask_, 4*asadIdx+agetIdx, 1, (isChipActive?1u:0));
					}
				}
			}
			catch (const CCfg::CConfig::NotFound &)
			{}
		}
	}

	LOG_DEBUG() << "Mask of active AsAd boards: " << std::showbase << std::hex << asadMask_ << std::dec;
	LOG_DEBUG() << "Mask of active AGET chips: " << std::showbase << std::hex << agetMask_ << std::dec;
}
//__________________________________________________________________________________________________
/**
 * Selects a unique AsAd board for a slow control command.
 * @param asadIdx Index of the AsAd board to select.
 */
void TestController::selectAsAdForSlowControl(const size_t & asadIdx)
{
	const uint16_t mask = (1ul << asadIdx);
	eccClient().setAsAdSlowControlMask(mask);
	eccClient().selectAsAdForRead(asadIdx);
}
//__________________________________________________________________________________________________
/**
 * Selects a unique AGET chip for a slow control command.
 * @param asadIdx Index of the AsAd board.
 * @param asadIdx Index of the AGET chip to select.
 */
void TestController::selectAgetForSlowControl(const size_t & asadIdx, const size_t agetIdx)
{
	const uint16_t mask = (1ul << (4*asadIdx+agetIdx));
	eccClient().setAgetSlowControlMask(mask);
	eccClient().selectAgetForRead(asadIdx, agetIdx);
}
//__________________________________________________________________________________________________
/**
 * Establishes connection to ECC server by allocating EccClient if not yet done.
 * @return Reference to allocated EccClient object.
 */
EccClient& TestController::eccClient()
{
	if (!eccClient_.get())
	{
		eccClient_.reset(new EccClient(testSummary().eccEndpoint().toString()));
		eccClient_->setVerbose(isVerbose());
		eccClient_->setActiveAsAdMask(asadMask_);
		eccClient_->setActiveAGetMask(agetMask_);
		LOG_INFO() << "Created client for connection to ECC server " << testSummary().eccEndpoint().toString();
		MSG_IF(verbose) << "Verbose mode is " << (isVerbose()?"":"NOT ") << "enabled";
		MSG_IF(verbose) << "Data display is " << (ControllerOptions::instance().isDataDisplayEnabled()?"":"NOT ") << "enabled";
	}
	return *eccClient_;
}
//__________________________________________________________________________________________________
mt::AlarmCallbackPrx TestController::alarmLogger()
{
	const ::utl::net::SocketAddress loggerEndpoint = ControllerOptions::instance().alarmLoggerEndpoint();
	Ice::ObjectAdapterPtr adapter = ::mdaq::utl::Server::find(loggerEndpoint).adapter();
	Ice::Identity id = adapter->getCommunicator()->stringToIdentity("AlarmLogger");

	Ice::ObjectPtr object = adapter->find(id);
	mt::AlarmCallbackPrx logger = 0;
	if (object)
	{
		logger = mt::AlarmCallbackPrx::checkedCast(adapter->createProxy(id));
	}
	return logger;
}
//__________________________________________________________________________________________________
void TestController::subscribeToAlarms(const size_t & coboId)
{
	try
	{
		// Check whether to monitor AsAd alarms
		bool checkAsAd = true;
		try
		{
			checkAsAd = coboConfig(coboId)("Alarms")("checkAsAd").getBoolValue();
		}
		catch (const CCfg::CConfig::NotFound & e)
		{
			LOG_WARN() << e.what();
		}
		// Enable monitoring of AsAd alarms
		LOG_INFO() << "Enabling periodic check of AsAd alarm bits by CoBo";
		eccClient().ecc()->setAsAdAlarmMonitoringEnabled(checkAsAd);
		// Remove all subscribers
		mt::AlarmServicePrx alarmSvce = eccClient().ecc()->alarmService();
		alarmSvce->reset();
		// Check whether to subscribe
		bool subscribe = true;
		try
		{
			subscribe = coboConfig(coboId)("Alarms")("subscribe").getBoolValue();
		}
		catch (const CCfg::CConfig::NotFound & e)
		{
			LOG_WARN() << e.what();
		}
		// Subscribe
		if (subscribe)
		{
			LOG_INFO() << "Subscribing to CoBo alarm notifications";
			const std::string id = std::string("GetController-") + QUuid::createUuid().toString().toStdString();
			alarmSvce->subscribe(id, alarmLogger());
		}
	}
	catch (::mdaq::utl::CmdException & e)
	{
		LOG_WARN() << e.reason;
	}
}
//__________________________________________________________________________________________________
bool TestController::isRunFinished()
{
	QMutexLocker guard(&mutex);
	return isRunFinished_;
}
//__________________________________________________________________________________________________
void TestController::setRunFinished(bool enabled)
{
	QMutexLocker guard(&mutex);
	isRunFinished_ = enabled;
}
//__________________________________________________________________________________________________
/**
 * Creates output files.
 */
void TestController::createOutputFiles(const std::string & prefix)
{
	// Create sub-directory
	QString outputPath = QString::fromStdString(testSummary().dataPath());
	if (not QDir(outputPath).exists())
	{
		LOG_INFO() << "Creating directory '" << outputPath.toStdString() << '\'';
		QDir().mkpath(outputPath);
	}

	// Create file for raw frames
	isSaveRawFrameEnabled_ = ControllerOptions::instance().isSaveRawFrameEnabled();
	if (isSaveRawFrameEnabled_)
	{
		createRawFrameFile(prefix, true);
	}
}
//_________________________________________________________________________________________________
void TestController::createTimeStamp()
{
	fileTimestamp = get::utl::buildTimeStamp();
}
//_________________________________________________________________________________________________
void TestController::createRawFrameFile(const std::string & prefix, bool newTimeStamp)
{
	closeRawFrameFile();

	// File timestamp and index
	if (newTimeStamp or fileTimestamp.empty())
	{
		createTimeStamp();
		rawFileIndex = 0;
	}
	else
	{
		rawFileIndex++;
	}
	rawFrameIndex = 0;
	// File prefix
	if (not prefix.empty())
	{
		filePrefix = prefix;
	}
	// Filename
	std::ostringstream oss;
	oss << testSummary().dataPath() << "/" << filePrefix << '_' << fileTimestamp << '_' << std::setfill('0') << std::setw(4) << rawFileIndex << ".graw";
	rawFrameFileName = oss.str();

	LOG_INFO() << "Creating new file: '" << rawFrameFileName << '\'';
	rawFrameFile.open(rawFrameFileName.c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
	if (not rawFrameFile)
	{
		LOG_ERROR() << "Could not create file '" << rawFrameFileName << "'.";
	}
	closeRawFrameFile();
}
//__________________________________________________________________________________________________
void TestController::closeRawFrameFile()
{
	if (rawFrameFile.is_open())
		rawFrameFile.close();
}
//__________________________________________________________________________________________________
void TestController::openRawFrameFile()
{
	if (not rawFrameFile.is_open())
	{
		rawFrameFile.open(rawFrameFileName.c_str(), std::ios::out | std::ios::app);
		if (rawFrameFile.is_open())
		{
			LOG_DEBUG() << "Appending data in file: " << rawFrameFileName;
		}
		else
		{
			LOG_ERROR() << "Could not open file: " << rawFrameFileName;
		}
	}
}
//__________________________________________________________________________________________________
void TestController::writeHeader(mfm::Frame & frame, std::ostream & stream)
{
	stream << "Header:\n";
	stream << std::dec;

	stream << "blockSize=" << frame.header().blkSize_B();
	stream << ' ';
	stream << "frameSize=" << frame.header().frameSize_B();
	stream << ' ';
	stream << "frameType=" << frame.findFormat().formatName();
	stream << ' ';
	stream << "revision=" << (short) frame.header().revision();
	stream << ' ';
	stream << "headerSize=" << frame.header().headerSize_B();
	stream << std::endl;

	stream << "itemSize=" << frame.header().itemSize_B();
	stream << ' ';
	stream << "itemCount=" << frame.header().itemCount();
	stream << ' ';
	stream << "eventTime=" << frame.headerField("eventTime").value< uint64_t >();
	stream << ' ';
	stream << "eventIdx=" << frame.headerField("eventIdx").value< uint32_t >();
	stream << std::endl;

	stream << "coboIdx=" << (short) frame.headerField("coboIdx").value< uint8_t >();
	stream << ' ';
	stream << "asadIdx=" << (short) frame.headerField("asadIdx").value< uint8_t >();
	stream << ' ';
	stream << "readOffset=" << frame.headerField("readOffset").value< uint16_t >();
	stream << std::endl;

	for (std::size_t agetIdx=0; agetIdx < 4; agetIdx++)
	{
		std::ostringstream oss;
		oss << "hitPat_" << agetIdx;
		mfm::DynamicBitset pattern;
		frame.headerField(oss.str()).value(pattern);
		stream << oss.str() << ": " << pattern;
		stream << std::endl;
		oss.str("");
		oss << "multip_" << agetIdx;
		stream << oss.str() << ": " << frame.headerField(oss.str()).value< uint16_t >();
		stream << std::endl;
	}
	// Print sliding window
	stream << "windowOut: " << frame.headerField("windowOut").value< uint32_t > () << std::endl;
	// Print last cell index
	try
	{
		for (short agetIdx=0; agetIdx < 4; ++agetIdx)
		{
			std::ostringstream oss;
			oss << "lastCell_" << agetIdx;
			stream << "lastCell_" << agetIdx << ": " << frame.headerField(oss.str()).value< uint16_t >() << std::endl;
		}
	}
	catch (const mfm::FieldNotFound &)
	{}
}
//__________________________________________________________________________________________________
/**
 * Writes the raw frame into the dedicated file.
 */
void TestController::writeFrameData(mfm::Frame & frame)
{
	// Store frame in file
	openRawFrameFile();
	frame.write(rawFrameFile);
	//closeRawFrameFile();
	rawFrameIndex++;
}
//__________________________________________________________________________________________________
/**
 *  calcul du seuil 4 bit pour un registre 64 bits
 */
uint64_t TestController::setChannelThreshold(uint64_t value, uint32_t n)
{
	uint64_t val = 0;
	for (uint32_t i = 0; i < n; i++)
	{
		val = val + (value << (i * 4));
	}
	LOG_INFO() << "Threshold value: " << QString("%1").arg(val, 64, 2).toStdString();
	return val;
}
//__________________________________________________________________________________________________
/**
 *  Sets the 4-bit individual channel thresholds in AGET registers 'reg8' and 'reg9'.
 */
void TestController::setChanRegThreshold(const std::size_t & asadIdx, const std::size_t & agetIdx, const uint64_t & value)
{
	uint64_t value1 = setChannelThreshold(value, 16);
	eccClient().agetRegWSC(asadIdx, agetIdx, "reg8", value1);
	eccClient().agetRegWSC(asadIdx, agetIdx, "reg8msb", value1);
	eccClient().agetRegWSC(asadIdx, agetIdx, "reg9", value1);
	eccClient().agetRegWSC(asadIdx, agetIdx, "reg9msb", value1);
}
//__________________________________________________________________________________________________
/**
 * Adjusts the 3-bit global threshold in AGET register 1.
 */
void TestController::setGlobalRegThreshold(const std::size_t & asadIdx, const std::size_t & agetIdx, const uint32_t & value)
{
	MSG_IF(verbose) << "Setting global threshold for AGET no. " << asadIdx*4+agetIdx << " to " << value;

	uint64_t regValue = eccClient().agetRegRSC(asadIdx, agetIdx, "reg1");
	::utl::BitFieldHelper< uint64_t >::setField(regValue, 19, 3, value);

	eccClient().agetRegWSC(asadIdx, agetIdx, "reg1", regValue);
}
//__________________________________________________________________________________________________
/**
 * Sets the 4 bits of AGET register 1 (state1<3:6>) corresponding to the shaper peaking time value.
 */
void TestController::setPeakingTimeIndex(const uint32_t & value)
{
	const uint32_t peakingTimeIndex = value % 16;
	LOG_INFO() << "Selecting shaper peaking time value no. " << peakingTimeIndex;

	// Find first active chip of the 16 chips
	size_t refChip;
	for (refChip=0; refChip < 16; ++refChip)
	{
		if (BitFieldHelper< uint16_t >::getBit(agetMask_, refChip))
			break;
	}
	const size_t agetIdx = refChip%4;
	const size_t asadIdx = refChip/4;
	uint32_t regValue = eccClient().agetRegRSC(asadIdx, agetIdx, "reg1");

	BitFieldHelper<uint32_t>::setField(regValue, 3, 4, peakingTimeIndex);
	eccClient().agetRegWSC("reg1", regValue);
}
//__________________________________________________________________________________________________
/**
 * Sets the 4 bits of AGET register 1 (state1<3:6>) corresponding to the shaper peaking time value.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 */
void TestController::setPeakingTimeIndex(const std::size_t & asadIdx, const std::size_t & agetIdx, const uint32_t & value)
{
	const uint32_t peakingTimeIndex = value % 16;
	MSG_IF(verbose) << "Selecting shaper peaking time value no. " << peakingTimeIndex;
	uint32_t regValue = eccClient().agetRegRSC(asadIdx, agetIdx, "reg1");

	BitFieldHelper<uint32_t>::setField(regValue, 3, 4, peakingTimeIndex);
	eccClient().agetRegWSC(asadIdx, agetIdx, "reg1", regValue);
}
//_________________________________________________________________________________________________
/**
 * Initializes test run loop.
 */
void TestController::initStopLoop()
{
	stopLoop = false;
}
//_________________________________________________________________________________________________
/**
 * Initialization of ECC server connection.
 * Initialization of the system for test algorithms
 * Creates data server and connects to target
 * Reduced-GET system parameters initialization.
 */
void TestController::loadAndPrepare()
{
	// Load hardware description
	eccClient().loadHwDescription(
			testSummary().targetEndpoint().toString(),
			testSummary().hwConfigPath());

	// Set 'initDone' bit to 'not configured' (0) state
	eccClient().setCoBoInitDone(false);
	// Tell CoBo board it will not be driven by Mutant
	eccClient().setCoBoMutantMode(false);

	// Check hardware type
	coboVersion_ = eccClient().checkHardwareVersion();

	// Check firmware version
	eccClient().checkFirmwareRelease();

	// Initialize firmware pipeline
	eccClient().pipeInit();

	// Set function of CoBo LEMO connectors
	setCoBoLemoModes(0);

	// Initialize trigger mode
	eccClient().triggerInit();

	// Initialize circular buffer
	eccClient().circularBufferInit();
	// Initialize AsAd registers
	asadInit(10);
	// Initialize CoBo and AsAd clock frequencies
	setCoBoWritingClockFrequency(0);
	asadSamplingClockInit(0);
	setCoBoReadingClockFrequency(0);
	asadReadingClockInit(0);
	// Set manual SCR data delay
	setReadDataDelay(0);
	// Set AsAd AGET input protections
	asadAgetInputManagerInit(0);
	// Initialize AGET registers
	agetInit();
	// Initialize AsAd ADC parameters
	eccClient().asadAdcInit(EccClient::DDR);
	// Initialize AsAd DAC registers
	eccClient().asadDacInit();
	// Initialize AsAd auto-protection
	asadAlertsInit(0);
	// Subscribe to alarm notifications
	subscribeToAlarms(0);
	asadMonitoringRead();
	// Initialize value of AGET global threshold (why now?)
	initializeGlobalThresholds(0, 2);
	// Compute automatic SCR data delay
	calibrateAutoReadDataDelay(0);
	// Select automatic or manual SCR data delay
	selectReadDataDelay(0);

	eccClient().sleep_ms(500);
}
//__________________________________________________________________________________________________
void TestController::unprepare()
{
	try
	{
		// Reconnect to node in case another client of the ECC changed current node
		eccClient().connectNode(testSummary().targetEndpoint().toString());

		// Power off AsAd boards
		eccClient().asadPowerOff();

		// Disconnect CoBo from data router
		daqDisconnect();
	}
	catch (::mdaq::utl::CmdException &)
	{
	}
}
//__________________________________________________________________________________________________
/**
 * Creates connection from (Reduced-CoBo) board to GetController data router.
 */
void TestController::daqConnect()
{
	std::string senderFlowType(ControllerOptions::instance().dataFlowType());
	ba::to_lower(senderFlowType);
	std::string flowEndpoint = ControllerOptions::instance().dataFlowEndpoint().toString();

	eccClient().daqConnect(flowEndpoint, senderFlowType);
}
//__________________________________________________________________________________________________
void TestController::daqDisconnect()
{
	try
	{
		eccClient().daqDisconnect();
	}
	catch (const std::exception &)
	{
	}
}
//__________________________________________________________________________________________________
/**
 * Starts or stops transfer of data present in the embedded circular buffer.
 */
void TestController::daqStart(bool status)
{
	const int coboIdx = 0; // FIXME: adapt to multiCoBo
	bool resetTime = true;
	try
	{
		resetTime = not coboConfig(coboIdx)("Module")("useExternalTimeReset").getBoolValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		// Keep default value
	}
	// Reconnect to node in case another client of the ECC changed current node
	eccClient().connectNode(testSummary().targetEndpoint().toString());

	eccClient().daqCtrl(status, resetTime);
}
//__________________________________________________________________________________________________
/**
 * Triggers automatic calibration of data delays, if required by configuration.
 */
void TestController::calibrateAutoReadDataDelay(int coboIdx)
{
	if (isCoboActive(coboIdx))
	{
		bool calcEnabled = true;
		try
		{
			calcEnabled = coboConfig(coboIdx)("Module")("calibrateAutoReadDataDelay").getBoolValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			// Keep default value
		}
		if (calcEnabled)
			eccClient().calibrateAutoReadDataDelay();
	}
}
//__________________________________________________________________________________________________
/**
 * Initializes the 3-bit global threshold of all active AGET chips. (Why?)
 */
void TestController::initializeGlobalThresholds(const size_t & coboIdx, const uint16_t value)
{
	LOG_INFO() << "Initializing global thresholds";

	eccClient().setAgetSlowControlMask(agetMask_);

	for (size_t asadIdx=0; asadIdx < 4; ++asadIdx)
	{
		for (size_t agetIdx=0; agetIdx < 4; ++agetIdx)
		{
			if (not isAgetActive(coboIdx, asadIdx, agetIdx)) continue;

			setGlobalRegThreshold(asadIdx, agetIdx, value);
		}
	}
}
//__________________________________________________________________________________________________
/**
 * Configures the test mode where the CoBo frames contain an incremental counter instead of data.
 */
void TestController::setAcqMode(const size_t & coboIdx)
{
	if (isCoboActive(coboIdx))
	{
		eccClient().acqTestMode(coboConfig(coboIdx)("Module")("isAcqTestMode").getBoolValue());
	}
}
//__________________________________________________________________________________________________
/**
 * Checks whether to use the automatically calibrated read data delays.
 *
 */
bool TestController::isSCRAutoDelayEnabled(const size_t & coboId, const size_t & asadId)
{
	bool enable = false;

	try
	{
		enable = asadConfig(coboId, asadId)("Control")("enableAutoReadDataDelay").getBoolValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		try
		{
			enable = coboConfig(coboId)("Module")("enableAutoReadDataDelay").getBoolValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			// Keep default value
		}
	}

	return enable;
}
//_________________________________________________________________________________________________
/**
 * Selects automatic or manual SCR data delay, depending on the value of the
 * Setup[Conditions].CoBo[*].Module.enableAutoReadDataDelay parameter
 * in the the configuration.
 * @param coboId Index of the CoBo board.
 */
void TestController::selectReadDataDelay(const size_t & coboId)
{
	if (isCoboActive(coboId))
	{
		for (size_t asadId = 0; asadId < 4u; asadId++ )
		{
			if (isAsadActive(coboId, asadId))
			{
				eccClient().asadSetAutoReadDataDelayEnabled(asadId, isSCRAutoDelayEnabled(coboId, asadId));
			}

		}
	}
}
//_________________________________________________________________________________________________
/**
 * Returns ID for channel chosen for single channel monitoring.
 * Used to set the fist 3 bits of the 'configuration2' AsAd register.
 */
uint8_t TestController::asadGetSingleChannelItem(const size_t & coboIdx, const size_t & asadIdx)
{
	std::string itemString = asadConfig(coboIdx, asadIdx)("Monitoring")("singleChannelItem").getStringValue("VDD1");
	// 000: VDD1; 001: Internal temperature; 010: External temperature; 100: VDD; 101: IDD
	uint8_t itemId = 0;
	if ("Tint" == itemString)
		itemId = 1;
	else if ("Text" == itemString)
		itemId = 2;
	else if ("VDD" == itemString)
		itemId = 4;
	else if ("IDD" == itemString)
		itemId = 5;
	return itemId;
}
//__________________________________________________________________________________________________
/**
 * Sets AsAd board auto-protection level.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring
 */
void TestController::asadSetAutoProtectionLevel(const size_t & coboIdx, const size_t & asadIdx)
{
	const bool enableAlarms = asadConfig(coboIdx, asadIdx)("Monitoring")("isAlarmSet").getBoolValue();

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
			flag = asadConfig(coboIdx, asadIdx)("Monitoring")("LimitsText")("disableAlertsTooHigh").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 2u, flag);
			flag = asadConfig(coboIdx, asadIdx)("Monitoring")("LimitsText")("disableAlertsTooLow").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 3u, flag);
			flag = asadConfig(coboIdx, asadIdx)("Monitoring")("LimitsText")("disableAlertsOutOfOrder").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 4u, flag);
			flag = asadConfig(coboIdx, asadIdx)("Monitoring")("LimitsVDD")("disableAlerts").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 6u, flag);
			flag = asadConfig(coboIdx, asadIdx)("Monitoring")("LimitsIDD")("disableAlerts").getBoolValue();
			::utl::BitFieldHelper< uint8_t >::setBit(interruptMask1, 7u, flag);
		}
		catch (const CCfg::CConfig::NotFound & e)
		{
			LOG_WARN() << e.what();
		}
	}

	// Leaving "interruptMask2" to its default value (has to do with VDD1)
	const uint8_t interruptMask2 = 0;

	eccClient().asadSetAutoProtectionLevel(interruptMask1, interruptMask2);
}
//__________________________________________________________________________________________________
/**
 * Sets AsAd monitoring device operational mode.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring
 */
void TestController::asadMonitoringSetOperationalMode(bool enableAlarms, bool averagingDisabled, bool singleChannelMode, uint8_t singleChannelItem, bool fastMode)
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

	eccClient().asadSetMonitoringMode(configuration1, configuration2, configuration3);

	LOG_DEBUG() << "Set operational mode of AsAd monitoring device to:" << std::hex << std::showbase
			<< " cfg1=" << (uint16_t) configuration1
			<< " cfg2=" << (uint16_t) configuration2
			<< " cfg3=" << (uint16_t) configuration3 << std::dec;
}
//_________________________________________________________________________________________________
/**
 * Defines AsAd auto-protection level and monitoring mode.
 */
void TestController::asadMonitoringSetOperationalMode(const size_t & coboIdx, const size_t & asadIdx)
{
	// Define operational mode
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
	bool fastMode = asadConfig(coboIdx, asadIdx)("Monitoring")("fastConversionEnabled").getBoolValue(false);

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
	bool singleChannelMode = asadConfig(coboIdx, asadIdx)("Monitoring")("singleChannelModeEnabled").getBoolValue(false);
	bool averagingDisabled = asadConfig(coboIdx, asadIdx)("Monitoring")("averagingDisabled").getBoolValue(false);
	uint8_t singleChannelItem = asadGetSingleChannelItem(coboIdx, asadIdx);

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
	const bool enableAlarms = asadConfig(coboIdx, asadIdx)("Monitoring")("isAlarmSet").getBoolValue();

	asadMonitoringSetOperationalMode(enableAlarms, fastMode, singleChannelItem, singleChannelMode, averagingDisabled);
}
//_________________________________________________________________________________________________
void TestController::asadMonitoringCheckVDD(const size_t & coboIdx, const size_t maxReads, const double delay_ms)
{
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
	{
		if (isAsadActive(coboIdx, asadIdx))
		{
			bool ok = false;
			const double vddMin = asadConfig(coboIdx, asadIdx)("Monitoring")("LimitsVDD")("lowerValue").getRealValue(3.2);
			LOG_INFO() << "Waiting for AsAd voltage in AsAd no. " << asadIdx << " to settle above " << vddMin << " V ...";
			for (size_t i=0; i < maxReads; ++i)
			{
				// Read VDD until high enough
				const double vdd = readAsAdVdd(asadIdx);
				if (vdd >= vddMin)
				{
					ok = true;
					break;
				}
				msleep(delay_ms);
			}
			if (not ok)
			{
				std::ostringstream msg;
				msg << "VDD voltage of AsAd board no. " << asadIdx << " did not reach " << vddMin
						<< " V! Check the hardware setup of your power supply or contact CENBG for further information.";
				throw mdaq::utl::CmdException(msg.str());
			}

		}
	}
}
//_________________________________________________________________________________________________
void TestController::asadMonitoringSetLimits(const size_t & coboIdx, const size_t & asadIdx)
{
	DBG_IF(verbose) << "Defining AsAd no. " << asadIdx << " monitoring limits";
	eccClient().asadMonitoringStop();

	asadSetExtTempLimits(coboIdx, asadIdx);
	asadSetVddLimits(coboIdx, asadIdx);
	asadSetIddLimits(coboIdx, asadIdx);
	//asadSetTintLimits(coboIdx, asadIdx);
	//asadSetVdd1Limits(coboIdx, asadIdx);
}
//_________________________________________________________________________________________________
double TestController::getVddLowLimit(int iCobo, int iAsad)
{
	double lowerVoltage = 0.0;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			lowerVoltage = asadConfig(iCobo, iAsad)("Monitoring")("LimitsVDD")("lowerValue").getRealValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			lowerVoltage = asadConfig(iCobo, iAsad)("Monitoring")("VddLowLimit").getRealValue();
		}
	}
	return lowerVoltage;
}
//_________________________________________________________________________________________________
double TestController::getVddHighLimit(int iCobo, int iAsad)
{
	double upperVoltage = 3.74;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			upperVoltage = asadConfig(iCobo, iAsad)("Monitoring")("LimitsVDD")("upperValue").getRealValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			upperVoltage = asadConfig(iCobo, iAsad)("Monitoring")("VddHighLimit").getRealValue();
		}
	}
	return upperVoltage;
}
//_________________________________________________________________________________________________
void TestController::asadSetVddLimits(const size_t & coboIdx, const size_t & asadIdx)
{
	const double lowerVoltage = getVddLowLimit(coboIdx, asadIdx);
	const double upperVoltage = getVddHighLimit(coboIdx, asadIdx);

	const uint8_t lowerLimit = 512. * lowerVoltage / 7.524;
	const uint8_t upperLimit = 512. * upperVoltage / 7.524;

	eccClient().wVddLimit(lowerLimit, upperLimit);
}
//_________________________________________________________________________________________________
double TestController::getVdd1LowLimit(int iCobo, int iAsad)
{
	double lowerVoltage = 0.0;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			lowerVoltage = asadConfig(iCobo, iAsad)("Monitoring")("LimitsVDD1")("lowerValue").getRealValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			lowerVoltage = asadConfig(iCobo, iAsad)("Monitoring")("Vdd1LowLimit").getRealValue();
		}
	}
	return lowerVoltage;
}
//_________________________________________________________________________________________________
double TestController::getVdd1HighLimit(int iCobo, int iAsad)
{
	double upperVoltage = 7.0;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			upperVoltage = asadConfig(iCobo, iAsad)("Monitoring")("LimitsVDD1")("upperValue").getRealValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			upperVoltage = asadConfig(iCobo, iAsad)("Monitoring")("Vdd1HighLimit").getRealValue();
		}
	}
	return upperVoltage;
}
//_________________________________________________________________________________________________
void TestController::asadSetVdd1Limits(const size_t & coboIdx, const size_t & asadIdx)
{
	const double lowerVoltage = getVdd1LowLimit(coboIdx, asadIdx);
	const double upperVoltage = getVdd1HighLimit(coboIdx, asadIdx);

	const uint8_t lowerLimit = 256. * (lowerVoltage - 0.06) / 7.;
	const uint8_t upperLimit = 256. * (upperVoltage - 0.06) / 7.;

	eccClient().wVdd1Limit(lowerLimit, upperLimit);
}
//_________________________________________________________________________________________________
double TestController::getIddLowLimit(int iCobo, int iAsad)
{
	double lowerCurrent = 0.0;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			lowerCurrent = asadConfig(iCobo, iAsad)("Monitoring")("LimitsIDD")("lowerValue").getRealValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			lowerCurrent = asadConfig(iCobo, iAsad)("Monitoring")("IddLowLimit").getRealValue();
		}
	}
	return lowerCurrent;
}
//_________________________________________________________________________________________________
double TestController::getIddHighLimit(int iCobo, int iAsad)
{
	double upperCurrent = 3.0;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			upperCurrent = asadConfig(iCobo, iAsad)("Monitoring")("LimitsIDD")("upperValue").getRealValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			upperCurrent = asadConfig(iCobo, iAsad)("Monitoring")("IddHighLimit").getRealValue();
		}
	}
	return upperCurrent;
}
//_________________________________________________________________________________________________
void TestController::asadSetIddLimits(const size_t & coboIdx, const size_t & asadIdx)
{
	const double lowerCurrent = getIddLowLimit(coboIdx, asadIdx);
	const double upperCurrent = getIddHighLimit(coboIdx, asadIdx);

	const uint8_t lowerLimit = 193.28 * lowerCurrent / 2.28;
	const uint8_t upperLimit = 193.28 * upperCurrent / 2.28;

	eccClient().wIddLimit(lowerLimit, upperLimit);
}
//_________________________________________________________________________________________________
uint8_t TestController::getIntTempLowLimit(int iCobo, int iAsad)
{
	uint8_t lowerValue = 0;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			lowerValue = asadConfig(iCobo, iAsad)("Monitoring")("LimitsTint")("lowerValue").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			lowerValue = asadConfig(iCobo, iAsad)("Monitoring")("InternalTempLowLimit").getIntValue();
		}
	}
	return lowerValue;
}
//_________________________________________________________________________________________________
uint8_t TestController::getIntTempHighLimit(int iCobo, int iAsad)
{
	uint8_t upperValue = 0xFF;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			upperValue = asadConfig(iCobo, iAsad)("Monitoring")("LimitsTint")("upperValue").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			upperValue = asadConfig(iCobo, iAsad)("Monitoring")("InternalTempHighLimit").getIntValue();
		}
	}
	return upperValue;
}
//_________________________________________________________________________________________________
void TestController::asadSetTintLimits(const size_t & coboIdx, const size_t & asadIdx)
{
	uint8_t lowerLimit = getIntTempLowLimit(coboIdx, asadIdx);
	uint8_t upperLimit = getIntTempHighLimit(coboIdx, asadIdx);

	eccClient().wIntTempLimit(lowerLimit, upperLimit);
}
//_________________________________________________________________________________________________
uint8_t TestController::getExtTempLowLimit(int iCobo, int iAsad)
{
	uint8_t lowerLimit = 0;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			lowerLimit = asadConfig(iCobo, iAsad)("Monitoring")("LimitsText")("lowerValue").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			lowerLimit = asadConfig(iCobo, iAsad)("Monitoring")("ExternalTempLowLimit").getIntValue();
		}
	}
	return lowerLimit;
}
//_________________________________________________________________________________________________
uint8_t TestController::getExtTempHighLimit(int iCobo, int iAsad)
{
	uint8_t upperLimit = 0xFF;
	if (isAsadActive(iCobo, iAsad))
	{
		try
		{
			upperLimit = asadConfig(iCobo, iAsad)("Monitoring")("LimitsText")("upperValue").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			upperLimit = asadConfig(iCobo, iAsad)("Monitoring")("ExternalTempHighLimit").getIntValue();
		}
	}
	return upperLimit;
}
//_________________________________________________________________________________________________
void TestController::asadSetExtTempLimits(const size_t & coboIdx, const size_t & asadIdx)
{
	uint8_t lowerLimit = getExtTempLowLimit(coboIdx, asadIdx);
	uint8_t upperLimit = getExtTempHighLimit(coboIdx, asadIdx);

	eccClient().wExtTempLimit(lowerLimit, upperLimit);
}
//__________________________________________________________________________________________________
/*
 * Checks in configuration file whether to set AsAd monitoring alarms ON.
 */
bool TestController::getMonitorAlarm(const size_t & coboIdx, const size_t & asadIdx)
{
	if (isAsadActive(coboIdx, asadIdx))
	{
		return (uint8_t) asadConfig(coboIdx, asadIdx)("Monitoring")("isAlarmSet").getBoolValue();
	}
	return false;
}
//__________________________________________________________________________________________________
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
void TestController::asadMonitoringInit(const size_t & coboIdx)
{
	LOG_DEBUG() << "Initializing AsAd monitoring registers";
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
	{
		if (isAsadActive(coboIdx, asadIdx))
		{
			// Select AsAd board for slow control
			selectAsAdForSlowControl(asadIdx);

			// Reset AsAd monitoring device
			eccClient().asadMonitoringReset();

			// Disable auto-protection
			eccClient().asadDisableAutoProtection(asadIdx);

			// Define upper and lower limits
			asadMonitoringSetLimits(coboIdx, asadIdx);

			// Define operational mode
			const bool fastMode = asadConfig(coboIdx, asadIdx)("Monitoring")("fastConversionEnabled").getBoolValue(false);
			const uint8_t singleChannelItem = asadGetSingleChannelItem(coboIdx, asadIdx);
			const bool singleChannelMode = asadConfig(coboIdx, asadIdx)("Monitoring")("singleChannelModeEnabled").getBoolValue(false);
			const bool averagingDisabled = asadConfig(coboIdx, asadIdx)("Monitoring")("averagingDisabled").getBoolValue(false);
			const bool enableAlarms = false;
			asadMonitoringSetOperationalMode(enableAlarms, fastMode, singleChannelItem, singleChannelMode, averagingDisabled);
		}
	}
	eccClient().setAsAdSlowControlMask(asadMask_);
}
//__________________________________________________________________________________________________
/**
 * Initializes AsAd monitoring registers.
 * See GET-AS-002-0004: AsAd Power Supply and Monitoring document.
 */
void TestController::asadAlertsInit(const size_t & coboIdx)
{
	LOG_DEBUG() << "Configuring AsAd auto-protection...";
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		if (isAsadActive(coboIdx, asadIdx))
		{
			// Select AsAd board for slow control
			selectAsAdForSlowControl(asadIdx);

			// Step 1: define upper and lower limits
			asadMonitoringSetLimits(coboIdx, asadIdx);
			// Step 2: define auto-protection level
			asadSetAutoProtectionLevel(coboIdx, asadIdx);
			// Step 3: define operational mode
			asadMonitoringSetOperationalMode(coboIdx, asadIdx);
		}
	}
	eccClient().setAsAdSlowControlMask(asadMask_);
}
//_________________________________________________________________________________________________
/**
 * Read AsAd voltages, currents, temperatures.
 */
void TestController::asadMonitoringRead(const size_t & coboIdx)
{
	eccClient().setAsAdSlowControlMask(asadMask_);
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
	{
		if (isAsadActive(coboIdx, asadIdx))
		{
			// Read AsAd voltages, currents, temperatures
			readAsAdVdd(asadIdx);
			readAsAdVdd1(asadIdx);
			readAsAdIdd(asadIdx);
			const int16_t Ti = readAsAdIntTemp(asadIdx);
			LOG_DEBUG() << "Read AsAd " << asadIdx << " temperature far from ADC: " << Ti << " deg C";
			const int16_t Te = readAsAdExtTemp(asadIdx);
			LOG_DEBUG() << "Read AsAd " << asadIdx << " temperature close to ADC: " << Te << " deg C";
		}
	}
}
//_________________________________________________________________________________________________
/**
 * Initialize AsAd registers
 */
void TestController::asadInit(const int & delay_ms, const size_t & coboIdx)
{
	// Check whether to abort if AsAd power supply connection is not detected
	uint8_t abortIfNoPowerSupply = 0xF;
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		BitFieldHelper< uint8_t >::setBit(abortIfNoPowerSupply, asadIdx, asadConfig(coboIdx, asadIdx)("Control")("checkPowerSupply").getBoolValue(true));
	}
	// Power on
	eccClient().asadInit(delay_ms, abortIfNoPowerSupply);

	eccClient().setAsAdSlowControlMask(asadMask_);

	// Calibrates the AsAd slow control
	eccClient().asadCal();

	// Configure and start AsAd monitoring (without auto-protection)
	asadMonitoringInit(0);

	// Check AsAd VDD voltages
	asadMonitoringCheckVDD(0);

	// Set AGET input protections
	eccClient().asadSetAgetInputProtections();
}
//_________________________________________________________________________________________________
/**
 * Reads AsAd VDD1 voltage (see AsAd Power Supply and Monitoring document).
 * Converts straight binary value to Volt.
 * @return Value of VDD1 voltage, in Volts.
 */
double TestController::readAsAdVdd1(const std::size_t & asadIdx)
{
	uint64_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		binValue = eccClient().readVdd1(asadIdx);
	}
	double value = (binValue != 0) ? (60e-3 + 7 * ((double) binValue) / 256) : 0.0;
	LOG_INFO() << "Read AsAd " << asadIdx << " VDD1 voltage: " << value << " V";
	return value;
}
//_________________________________________________________________________________________________
/**
 * Reads AsAd VDD voltage (see AsAd Power Supply and Monitoring document).
 * Converts straight binary value to Volt.
 * @return Value of VDD voltage, in Volts.
 */
double TestController::readAsAdVdd(const std::size_t & asadIdx)
{
	uint64_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		binValue = eccClient().readVdd(asadIdx);
	}
	double value = (binValue != 0) ? (2.28 / 256) * (3.3 / 2) * binValue : 0.0;
	LOG_INFO() << "Read AsAd " << asadIdx << " VDD voltage: " << value << " V";
	return value;
}
//_________________________________________________________________________________________________
/**
 * Reads AsAd IDD current (see AsAd Power Supply and Monitoring document).
 * Converts straight binary value to A.
 * @return Value of IDD current, in A.
 */
double TestController::readAsAdIdd(const std::size_t & asadIdx)
{
	uint64_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		binValue = eccClient().readIdd(asadIdx);
	}
	double value = binValue != 0 ? (2.28 / 256) * (1 / 0.755) * binValue : 0.0;
	LOG_INFO() << "Read AsAd " << asadIdx << " IDD current: " << value << " A";
	return value;
}
//_________________________________________________________________________________________________
/**
 * Reads AsAd on-chip temperature (see AsAd Power Supply and Monitoring document).
 * Converts straight binary value to degrees C.
 * @return Value of on-chip temperature, in degrees C.
 */
int8_t TestController::readAsAdIntTemp(const std::size_t & asadIdx)
{
	uint8_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		binValue = eccClient().readIntTemp(asadIdx);
	}
	uint8_t value = ((binValue & 0x80) == 0x8)? ~(binValue - 1) : binValue;
	return value;
}
//_________________________________________________________________________________________________
/**
 * Reads AsAd ADC-chip surroundings temperature (see AsAd Power Supply and Monitoring document).
 * Converts straight binary value to degrees C.
 * @return Value of ADC-chip surroundings temperature, in degrees C.
 */
int8_t TestController::readAsAdExtTemp(const std::size_t & asadIdx)
{
	uint8_t binValue = 0;
	for (int i = 0; i < 10 and binValue == 0; ++i)
	{
		binValue = eccClient().readExtTemp(asadIdx);
	}
	uint8_t value = ((binValue & 0x80) == 0x8) ? ~(binValue - 1) : binValue;
	return value;
}
//_________________________________________________________________________________________________
/**
 * Sets AsAd register for inspection manager to select signals.
 */
void TestController::setAsAdInspectionManager(const std::size_t & coboIdx, const std::size_t & asadIdx)
{
	std::string device, isp1, isp2;
	if (getInspectionLinesParameters(device, isp1, isp2, coboIdx, asadIdx))
	{
		LOG_INFO() << "Setting AsAd no. " << asadIdx << " inspection manager for signals (" << isp1 << ", " << isp2 << ") of device " << device;
		eccClient().asadSetISPM(coboIdx, asadIdx, device, isp1, isp2);
	}
}
//_________________________________________________________________________________________________
/**
 * Initializes CoBo registers to communicate with AGET
 */
void TestController::agetInit()
{
	eccClient().selectCalibChips(rc::CoBoConfigParser::getCalibrationChips(coboConfig(0)));
	eccClient().agetInit();
}
//_________________________________________________________________________________________________
/**
 * Starts a software triggered read.
 */
void TestController::triggerNow(int iCobo)
{
	if (iCobo == 0)
	{
		eccClient().triggerNow();
	}
}
//_________________________________________________________________________________________________
/**
 * Configures the activation of the AsAd ADC channels and the AsAd AGET Inputs Manager.
 * See Texas Instruments ADS6422 Data Sheet.
 * See AsAd Slow Control Protocol & Registers Mapping
 * See AsAd Switched Capacitor Array Management, "AGET Inputs Manager".
 *
 *  ADC register A (ADC0):
 * 	D10	<RST> 		S/W RESET
 * 	D9	0
 * 	D8	0
 * 	D7	0
 * 	D6	0
 * 	D5	<REF> 		INTERNAL OR EXTERNAL
 * 	D4	<PDN CHD>	POWER DOWN CH D
 * 	D3	<PDN CHC>	POWER DOWN CH C
 * 	D2	<PDN CHB>	POWER DOWN CH B
 * 	D1	<PDN CHA>	POWER DOWN CH A
 * 	D0	<PDN>		GLOBAL POWER DOWN
 *
 * Input Manager IM0 register:
 *  SWZ(3..0) Dbi(3..0)
 *  SWZ(x) = 1	Connect ground to ZAP
 *  SWZ(x) = 0	Connect VDD to ZAP
 *
 *  Dbi(x) = 1	Enable anti-sparkling protection for AGET
 *  Dbi(x) = 0	Disable anti-sparkling protection for AGET
 *
 * Input Manager IM1 register:
 *  ENW(3..0) ENR(3..0)
 *  ENW(x) = 1	Enable WCK clock at the input of AGET x
 *  			Enable SCA Write signal to AGETx
 *  ENW(x) = 0	Disable WCK clock at the input of AGET x
 *  			Disable SCA Write signal to AGETx
 *  			e.g. disable sampling command
 *  ENR(x) = 1	Enable RCK clock at the input of AGET x
 *  			Enable SCA Read signal to AGETx
 *  ENR(x) = 0	Disable RCK clock at the input of AGET x
 *  			Disable SCA Read signal to AGETx
 *  			e.g. disable readout command
 */
void TestController::asadAgetInputManagerInit(const size_t & coboIdx)
{
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		if (not isAsadActive(coboIdx, asadIdx)) continue;
		MSG_IF(verbose) << "Activating AsAd no. " << asadIdx << " ADC channels";

		// Default values
		uint8_t adc0 = 0;	// Enable all four ADC channels
		uint8_t im0 = 0xF; 	// Connect VDD to ZAP, Enable anti-sparkling protection
		uint8_t im1 = 0xFF;	// Allow sampling and readout for all four AGET chips

		// Loop over AGET chips
		for (size_t agetIdx=0; agetIdx < 4u; ++agetIdx)
		{
			// Check whether to activate slot
			bool isChipActive = isAgetActive(coboIdx, asadIdx, agetIdx);
			bool enablePosition = true;
			bool connectGroundToZAP = false; // SWZ
			bool antiSparklingEnabled = true; // Dbi

			try
			{
				enablePosition = agetConfig(coboIdx, asadIdx, agetIdx)("Control")("enableAsAdPosition").getBoolValue();
				connectGroundToZAP = agetConfig(coboIdx, asadIdx, agetIdx)("Control")("connectGroundToZAP").getBoolValue();
				antiSparklingEnabled = agetConfig(coboIdx, asadIdx, agetIdx)("Control")("antiSparklingEnabled").getBoolValue();
			}
			catch (const CCfg::CConfig::NotFound & e)
			{
				//LOG_WARN() << e.what();
			}

			// Disable slot
			if (not isChipActive or not enablePosition)
			{
				MSG_IF(verbose) << "Disabling position no. " << agetIdx << " of AsAd no. " << asadIdx;
				// Power down ADC channel
				BitFieldHelper< uint8_t >::setBit(adc0, 1 + agetIdx, true);
				// Disable sampling command on AGET input
				BitFieldHelper< uint8_t >::setBit(im1, 4 + agetIdx, false);
				// Disable readout command on AGET input
				BitFieldHelper< uint8_t >::setBit(im1, agetIdx, false);
			}

			// Connect/disconnect VDD to/from ZAP interface
			BitFieldHelper< uint8_t >::setBit(im0, 4 + agetIdx, connectGroundToZAP);
			// Enable/disable anti-sparkling protection for AGET
			BitFieldHelper< uint8_t >::setBit(im0, agetIdx, antiSparklingEnabled);
		}

		// Set mask for broadcast of AsAd slow control command
		selectAsAdForSlowControl(asadIdx);

		// Send command
		eccClient().asadConfigureAgetInputs(adc0, im0, im1);
	}

	eccClient().setAsAdSlowControlMask(asadMask_);
}
//_________________________________________________________________________________________________
uint16_t TestController::getInitGloablThresh(int iCobo, int iAsad, int iAget)
{
	uint16_t val = 0;
	if (isAgetActive(iCobo, iAsad, iAget))
	{
		val = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")(
				"GlobalThresholdValue").getIntValue();
		val = val * 16;
	}
	return val;
}
//_________________________________________________________________________________________________
/**
 * Analyzes the parameters in the configuration file and extracts values for the AGET 'reg1' register.
 */
uint64_t TestController::retrievalAget_reg1(int iCobo, int iAsad, int iAget, const uint16_t & agetVersion)
{
	uint64_t regValue = UINT64_C(0);

	if (isAgetActive(iCobo, iAsad, iAget))
	{
		MSG_IF(verbose) << "AGET selected  for 'reg1': " << agetConfig(iCobo, iAsad, iAget).getIndex();
		uint64_t fieldValue;

		// "Icsa". If set, the nominal CSA bias current is multiplied by two.
		const bool isIcsaSet = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("isIcsax2").getBoolValue();
		MSG_IF(verbose) << "isIcsax2: " << isIcsaSet;
		fieldValue = isIcsaSet?0x1:0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 0, 1, fieldValue);

		 // "Gain". LSb and MSb of the internal test capacitor in Test mode.
		const std::string testCapacitorValue = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("TestModeRange").getStringValue();
		MSG_IF(verbose) << "Internal test capacitor ('TestModeRange'): " << testCapacitorValue;
		// We expect Farad sub-units or, in old config files, Coulomb sub-units
		if (ba::starts_with(testCapacitorValue, "120f") or ba::starts_with(testCapacitorValue, "120 f"))
			fieldValue = UINT64_C(0x0);
		else if (ba::starts_with(testCapacitorValue, "240f") or ba::starts_with(testCapacitorValue, "240 f"))
			fieldValue = UINT64_C(0x1);
		else if (ba::starts_with(testCapacitorValue, "1p") or ba::starts_with(testCapacitorValue, "1 p") )
			fieldValue = UINT64_C(0x2);
		else if (ba::starts_with(testCapacitorValue, "10p") or ba::starts_with(testCapacitorValue, "10 p"))
			fieldValue = UINT64_C(0x3);
		else
		{
			LOG_WARN() << "Unexpected value for internal test capacitor: " << testCapacitorValue << ". Using 120 fF instead.";
			fieldValue = UINT64_C(0x0);
		}
		BitFieldHelper<uint64_t>::setField(regValue, 1, 2, fieldValue);

		// "Time". Sets the peaking time (ns) of the shaper by switching resistors on the PZC  SK filters.
		const int peakingTime = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("peackingTime").getIntValue();
		MSG_IF(verbose) << "Peaking time: " << peakingTime << " ns";
		fieldValue = rc::CoBoConfigParser::peakingTimeIndex(peakingTime, agetVersion);
		BitFieldHelper<uint64_t>::setField(regValue, 3, 4, fieldValue);

		// "Test". Defines the test mode: none (0), calibration (01), test (10), functionality (11).
		const std::string testMode = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("TestModeSelection").getStringValue();
		MSG_IF(verbose) << "TestModeSelection: " << testMode;
		if (testMode == "nothing")
			fieldValue = 0x0;
		else if (testMode == "calibration")
			fieldValue = 0x1;
		else if (testMode == "test")
			fieldValue = 0x2;
		else if (testMode == "functionality")
			fieldValue = 0x3;
		else
			fieldValue = 0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 7, 2, fieldValue);

		// "Integrator_mode". Bit integrator mode
		const std::string integratorMode = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("IntegrationMode").getStringValue();
		MSG_IF(verbose) << "IntegrationMode: " << integratorMode;
		if (integratorMode == "variable")
			fieldValue = 0x0;
		else
			fieldValue = 0x1;
		BitFieldHelper<uint64_t>::setField(regValue, 9, 1, fieldValue);

		// "SCApointer". Bit SCA pointer 0 to pointer 1.
		const std::string scaPointer = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("SCA_Pointer").getStringValue();
		MSG_IF(verbose) << "SCA_Pointer: " << scaPointer;
		if (scaPointer == "anySCAPointer")
			fieldValue = 0x3;
		else if (scaPointer == "col_0or128or256or384")
			fieldValue = 0x2;
		else if (scaPointer == "col_0or256")
			fieldValue = 0x1;
		else if (scaPointer == "column_0")
			fieldValue = 0x0;
		else
			fieldValue = 0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 10, 2, fieldValue);

		// "SCAsplitting". Bit SCA splitting; if set, enables the 2p mode by dividing the memory in two blocks.
		if (agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("SCA_Splitting").getBoolValue())
			fieldValue = 0x1;
		else
			fieldValue = 0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 12, 1, fieldValue);
		MSG_IF(verbose) << "SCA_Splitting: " << fieldValue;

		// "Mode32Channels". Bit 32 channels mode. If set, only the first 32 channels are used.
		if (agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("is32channels").getBoolValue())
			fieldValue = 0x1;
		else
			fieldValue = 0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 13, 1, fieldValue);
		MSG_IF(verbose) << "is32channels: " << fieldValue;

		// "Readout_mode". Bit readout mode; switches between full and partial readout modes. If set, full readout.
		if (agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("isAllChannelRead").getBoolValue())
			fieldValue = 0x1;
		else
			fieldValue = 0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 14, 1, fieldValue);
		MSG_IF(verbose) << "isAllChannelRead: " << fieldValue;

		// "FPNreadout". Bit FPN readout: if set, forces the readout of the FPN channels (e.g. even in partial readout mode).
		if (agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("isFPNRead").getBoolValue())
			fieldValue = 0x1;
		else
			fieldValue = 0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 15, 1, fieldValue);
		MSG_IF(verbose) << "isFPNRead: " << fieldValue;

		// "Polarity". Bit polarity. Controls the value of the DC voltage levels. Usually set to zero?
		if (agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("isPositivePolarity").getBoolValue())
			fieldValue = 0x1;
		else
			fieldValue = 0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 16, 1, fieldValue);
		MSG_IF(verbose) << "isPositivePolarity: " << fieldValue;

		// "Vicm". Input common mode voltage (VICM) of the analog output buffer. Usually 0 e.g. 1.25 V?
		const std::string vicmValue = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("Vicm").getStringValue();
		if (vicmValue == "1.25V")
			fieldValue = 0x0;
		else if (vicmValue == "1.35V")
			fieldValue = 0x1;
		else if (vicmValue == "1.55V")
			fieldValue = 0x2;
		else if (vicmValue == "1.65V")
			fieldValue = 0x3;
		else
			fieldValue = 0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 17, 2, fieldValue);
		MSG_IF(verbose) << "Vicm: " << vicmValue;

		// "DACthreshold".
		MSG_IF(verbose) << "GlobalThresholdValue";
		BitFieldHelper<uint64_t>::setField(regValue, 19, 3, agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")(
				"GlobalThresholdValue").getIntValue());

		// "DACpolarity".
		MSG_IF(verbose) << "isThresholdSignedPositive";
		fieldValue = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("isThresholdSignedPositive")
				.getBoolValue()?0x1:0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 22, 1, fieldValue);

		// "Trigger_veto". Trigger bit: specifies the use of the veto on the trigger building. Defines how long after a hit the channel can be read again.
		// Values: 0 (none), 01 ( 4 microseconds), 10 (hit register width).
		const std::string triggerVetoOption = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")(
				"TriggerVetoOption").getStringValue();
		if (triggerVetoOption == "none")
			fieldValue = 0x0;
		else if (triggerVetoOption == "4us")
			fieldValue = 0x1;
		else if (triggerVetoOption == "HitRegisterWidth")
			fieldValue = 0x2;
		else if (triggerVetoOption == "undefined")
			fieldValue = 0x3;
		else
			fieldValue = 0x0;
		MSG_IF(verbose) << "TriggerVetoOption: " << triggerVetoOption;
		BitFieldHelper<uint64_t>::setField(regValue, 23, 2, fieldValue);

		// "Synchro_discri".
		fieldValue = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("isckwriteSynchro")
				.getBoolValue()?0x1:0x0;
		MSG_IF(verbose) << "isckwriteSynchro: " << fieldValue;
		BitFieldHelper<uint64_t>::setField(regValue, 25, 1, fieldValue);

		// "tot".
		fieldValue = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("isTOTActive").
				getBoolValue()? 0x1:0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 26, 1, fieldValue);
		MSG_IF(verbose) << "isTOTActive: " << fieldValue;

		// "Range_trigg_width"
		fieldValue = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("isTriggerWidth200ns")
				.getBoolValue()?0x1:0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 27, 1, fieldValue);
		MSG_IF(verbose) << "isTriggerWidth200ns: " << fieldValue;

		const std::string triggerWidthRange = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")(
				"TriggerWidthRange").getStringValue();
		MSG_IF(verbose) << "TriggerWidthRange: " << triggerWidthRange;
		if (triggerWidthRange == "width1")
		{
			// "Lsb_trigg_width"
			BitFieldHelper<uint64_t>::setField(regValue, 28, 1, 0x1);
			// "Msb_trigg_width"
			BitFieldHelper<uint64_t>::setField(regValue, 29, 1, 0x0);
		}
		else if (triggerWidthRange == "width2")
		{
			// "Lsb_trigg_width"
			BitFieldHelper<uint64_t>::setField(regValue, 28, 1, 0x0);
			// "Msb_trigg_width"
			BitFieldHelper<uint64_t>::setField(regValue, 29, 1, 0x1);
		}
		else if (triggerWidthRange == "width3")
		{
			// "Lsb_trigg_width"
			BitFieldHelper<uint64_t>::setField(regValue, 28, 1, 0x1);
			// "Msb_trigg_width"
			BitFieldHelper<uint64_t>::setField(regValue, 29, 1, 0x1);
		}
		else // "width0"
		{
			// "Lsb_trigg_width"
			BitFieldHelper<uint64_t>::setField(regValue, 28, 1, 0x0);
			// "Msb_trigg_width"
			BitFieldHelper<uint64_t>::setField(regValue, 29, 1, 0x0);
		}

		// "External".
		const std::string externalLink = agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("ExternalLink").getStringValue();
		MSG_IF(verbose) << "ExternalLink: " << externalLink;
		if (externalLink == "none")
			fieldValue = 0x0;
		else if (externalLink == "SKFilterInput")
			fieldValue = 0x1;
		else if (externalLink == "Gain-2Input")
			fieldValue = 0x2;
		else if (externalLink == "CSAStandby")
			fieldValue = 0x3;
		else
			fieldValue = 0x0;
		BitFieldHelper<uint64_t>::setField(regValue, 30, 2, fieldValue);
	}
	return regValue;
}

/**
 * analyse le fichier de parametres et recupere les valeurs correspndants
 * au registre 2 de aget
 */
uint64_t TestController::retrievalAget_reg2(int iCobo, int iAsad, int iAget)
{
	uint64_t datareg = 0ull;
	if (isAgetActive(iCobo, iAsad, iAget))
	{
		MSG_IF(verbose) << "AGET selected  for 'reg2': " << agetConfig(iCobo,
					iAsad, iAget).getIndex();

		MSG_IF(verbose) << "DebugSelection";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("DebugSelection").getStringValue()
				== "standby"))
			datareg = datareg + (0x0);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"DebugSelection").getStringValue() == "CSA"))
			datareg = datareg + (0x1);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"DebugSelection").getStringValue() == "CR"))
			datareg = datareg + (0x2);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"DebugSelection").getStringValue() == "Gain-2"))
			datareg = datareg + (0x3);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"DebugSelection").getStringValue() == "PositiveInputDiscri"))
			datareg = datareg + (0x4);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"DebugSelection").getStringValue() == "NegativeInputDiscri"))
			datareg = datareg + (0x5);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"DebugSelection").getStringValue() == "trigger"))
			datareg = datareg + (0x6);
		else
			datareg = datareg + (0x7);

		MSG_IF(verbose) << "isReadFromColumn0";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isReadFromColumn0").getBoolValue()))
			datareg = datareg + (0x1 << 3);
		else
			datareg = datareg + (0x0 << 3);

		MSG_IF(verbose) << "isDigitalOutputTest";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isDigitalOutputTest").getBoolValue()))
			datareg = datareg + (0x1 << 4);
		else
			datareg = datareg + (0x0 << 4);

		MSG_IF(verbose) << "isResetLevelUndefined";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isResetLevelUndefined").getBoolValue()))
			datareg = datareg + (0x1 << 5);
		else
			datareg = datareg + (0x0 << 5);

		MSG_IF(verbose) << "isMarkerEnable";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("isMarkerEnable").getBoolValue()))
			datareg = datareg + (0x1 << 6);
		else
			datareg = datareg + (0x0 << 6);

		MSG_IF(verbose) << "isResetLevelVdd";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isResetLevelVdd").getBoolValue()))
			datareg = datareg + (0x1 << 7);
		else
			datareg = datareg + (0x0 << 7);

		MSG_IF(verbose) << "isBoostPower";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("isBoostPower").getBoolValue()))
			datareg = datareg + (0x1 << 8);
		else
			datareg = datareg + (0x0 << 8);

		MSG_IF(verbose) << "isOutputscCkSynchro";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isOutputscCkSynchro").getBoolValue()))
			datareg = datareg + (0x1 << 9);
		else
			datareg = datareg + (0x0 << 9);

		MSG_IF(verbose) << "isSynchroEdgeFalling";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isSynchroEdgeFalling").getBoolValue()))
			datareg = datareg + (0x1 << 10);
		else
			datareg = datareg + (0x0 << 10);

		MSG_IF(verbose) << "isSCoutputBufferInTriState";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isSCoutputBufferInTriState").getBoolValue()))
			datareg = datareg + (0x1 << 11);
		else
			datareg = datareg + (0x0 << 11);

		MSG_IF(verbose) << "currentBias";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("currentBias").getStringValue()
				== "current0"))
			datareg = datareg + (0x0 << 12);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"currentBias").getStringValue() == "current1"))
			datareg = datareg + (0x1 << 12);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"currentBias").getStringValue() == "current2"))
			datareg = datareg + (0x2 << 12);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"currentBias").getStringValue() == "current3"))
			datareg = datareg + (0x3 << 12);

		MSG_IF(verbose) << "bufferCurrent";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("bufferCurrent").getStringValue()
				== "1.503mA"))
			datareg = datareg + (0x0 << 14);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"bufferCurrent").getStringValue() == "1.914mA"))
			datareg = datareg + (0x1 << 14);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"bufferCurrent").getStringValue() == "2.7mA"))
			datareg = datareg + (0x2 << 14);
		else if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"bufferCurrent").getStringValue() == "4.870mA"))
			datareg = datareg + (0x3 << 14);

		MSG_IF(verbose) << "isPWDonWrite";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("isPWDonWrite").getBoolValue()))
			datareg = datareg + (0x1 << 16);
		else
			datareg = datareg + (0x0 << 16);

		MSG_IF(verbose) << "isPWDonRead";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("isPWDonRead").getBoolValue()))
			datareg = datareg + (0x1 << 17);
		else
			datareg = datareg + (0x0 << 17);

		MSG_IF(verbose) << "isAlternatePower";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isAlternatePower").getBoolValue()))
			datareg = datareg + (0x1 << 18);
		else
			datareg = datareg + (0x0 << 18);

		MSG_IF(verbose) << "is1bitShortResetReadSeq";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"is1bitShortResetReadSeq").getBoolValue()))
			datareg = datareg + (0x1 << 19);
		else
			datareg = datareg + (0x0 << 19);

		MSG_IF(verbose) << "isTriggerOutputDisable";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isTriggerOutputDisable").getBoolValue()))
			datareg = datareg + (0x1 << 20);
		else
			datareg = datareg + (0x0 << 20);

		MSG_IF(verbose) << "isRisingedgeWriteAutoresetBank";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"isRisingedgeWriteAutoresetBank").getBoolValue()))
			datareg = datareg + (0x1 << 21);
		else
			datareg = datareg + (0x0 << 21);

		MSG_IF(verbose) << "islvdsTriggerOutput";

		if ((agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")(
				"islvdsTriggerOutput").getBoolValue()))
			datareg = datareg + (0x1 << 22);
		else
			datareg = datareg + (0x0 << 22);

		// Bit "Input_dynamic_range"
		// The input dynamic range of the discriminator can be fixed to 5% or 17.5% of the input dynamic range of the analog channel.
		bool useInputDynamicRange5pct = false; // Default value of 0 for compatibility with AGET 2.1 where this bit does no exist.
		try
		{
			useInputDynamicRange5pct = agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("inputDynamicRange_5pct").getBoolValue();
		}
		catch (const CCfg::CConfig::NotFound&)
		{
			// keep default value
		}
		const uint64_t bitValue = useInputDynamicRange5pct ? 0x1 : 0x0;
		MSG_IF(verbose) << "Input_dynamic_range: " << bitValue;
		BitFieldHelper<uint64_t>::setField(datareg, 24, 1, bitValue);

		// Bits 25-26 "trigger_LVDS_level". Specific to ASTRE chip, added Nov 2016.
		// Controls multiplicity level of LVDS trigger, if bit 22 is set.
		uint64_t pow_lvds_level = (uint64_t) agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("triggerLvdsThreshold").getIntValue(1u);
		uint64_t lvds_level = logtwo(pow_lvds_level);
		MSG_IF(verbose) << "LVDS trigger threshold: " << lvds_level;
		BitFieldHelper<uint64_t>::setField(datareg, 25, 2, lvds_level);

		// Bits 27-28 "unity_multiplicity_level". Specific to ASTRE chip, added Nov 2016.
		// Controls unity multiplicity (e.g. per channel).
		// AGET behavior for 0; for 1, 32 channels cover all range; for 2, 16 channels; or 3, 8 channels.
		uint64_t pow_unit_mult = (uint64_t) agetConfig(iCobo, iAsad, iAget)("Global")("Reg2")("multiplicityLevel").getIntValue(1u);
		uint64_t unit_mult = logtwo(pow_unit_mult);
		MSG_IF(verbose) << "Multiplicity unity level: " << unit_mult;
		BitFieldHelper<uint64_t>::setField(datareg, 27, 2, unit_mult);

	}
	return datareg;
}

/**
 * @param reg Name of AGET register to fill.
 * @param iCobo CoBo index
 * @param iAsad AsAd index
 * @param iAget AGET index
 * @return Value of the specified register corresponding to configuration.
 */
uint64_t TestController::retrievalChan_reg(const std::string & reg, int iCobo,
		int iAsad, int iAget)
{
	uint64_t datareg = 0ull;
	if (isAgetActive(iCobo, iAsad, iAget))
	{
		MSG_IF(verbose) << "AGET selected for '" << reg << "': " << agetConfig(iCobo, iAsad, iAget).getIndex();

		CCfg::CConfig::Iterator itChan = agetConfig(iCobo, iAsad, iAget).iterate("channel");

		while (itChan.hasNext())
		{
			itChan.next();

			if (CCfg::DEFAULT_OBJ_INDEX == itChan->getIndex())
				continue;

			if ((*itChan)("isActive").getBoolValue())
			{
				uint64_t bitMask = 0ll;
				int iChan = boost::lexical_cast<int>(itChan->getIndex());

				if (reg == "reg3")
					bitMask = agetChanBitmask_reg3(iChan, itChan);
				else if (reg == "reg4")
					bitMask = agetChanBitmask_reg4(iChan, itChan);
				else if (reg == "reg6")
					bitMask = agetChanBitmask_reg6(iChan, itChan);
				else if (reg == "reg7")
					bitMask = agetChanBitmask_reg7(iChan, itChan);
				else if (reg == "reg8msb")
					bitMask = rc::CoBoConfigParser::agetChanBitmask("reg8msb", iChan, itChan);
				else if (reg == "reg8")
					bitMask = rc::CoBoConfigParser::agetChanBitmask("reg8", iChan, itChan);
				else if (reg == "reg9msb")
					bitMask = rc::CoBoConfigParser::agetChanBitmask("reg9msb", iChan, itChan);
				else if (reg == "reg9")
					bitMask = rc::CoBoConfigParser::agetChanBitmask("reg9", iChan, itChan);
				else if (reg == "reg10")
					bitMask = agetChanBitmask_reg10(iChan, itChan);
				else if (reg == "reg11")
					bitMask = agetChanBitmask_reg11(iChan, itChan);
				else
					throw std::string("Register ") + reg + " does not exist";

				datareg |= bitMask;

				//LOG_INFO() << "Channel: %1 --> 0x%2 --> 0x%3.").arg(iChan).arg(bitMask, 0, 16) .arg(datareg, 0, 16));
			}
		}
	}
	DBG_IF(verbose) << "Retrieved value '" << std::showbase << std::hex << datareg << std::dec << "' from config. for AGET register '" << reg << "'";
	return datareg;
}

/**
 * analyse le fichier de parametres et recupere les valeurs correspndants
 * au registre 3 de aget
 */
uint64_t TestController::agetChanBitmask_reg3(int iChan,
		CCfg::CConfig::Iterator& itChan)
{
	uint64_t bitMask = 0ull;
	if ((*itChan)("isSelectedforTestMode").getBoolValue())
	{
		if (iChan >= 0 and iChan <= 33)
		{
			int shiftIndex = 33 - iChan;
			bitMask = 1ull << shiftIndex;
		}
	}
	return bitMask;
}

/**
 * analyse le fichier de parametres et recupere les valeurs correspndants
 * au registre 4 de aget
 */
uint64_t TestController::agetChanBitmask_reg4(int iChan,
		CCfg::CConfig::Iterator& itChan)
{
	uint64_t bitMask = 0ull;
	if ((*itChan)("isSelectedforTestMode").getBoolValue())
	{
		if (iChan >= 34 and iChan <= 67)
		{
			int shiftIndex = iChan - 34;
			bitMask = 1ull << shiftIndex;
		}
	}
	return bitMask;
}

/**
 * analyse le fichier de parametres et recupere les valeurs correspndants
 * au registre 6 de aget
 */
uint64_t TestController::agetChanBitmask_reg6(int iChan,
		CCfg::CConfig::Iterator& itChan)
{
	uint64_t bitMask = 0ull;
	if (iChan >= 0 and iChan <= 33 and iChan != 11 and iChan != 22)
	{
		int shiftIndex = 62 - 2 * iChan + 4;
		if (iChan < 11)
		{
			shiftIndex = 62 - 2 * iChan;
		}
		else if (iChan < 22 and iChan > 11)
		{
			shiftIndex = 62 - 2 * iChan + 2;
		}

		if ((*itChan)("Gain").getStringValue() == "120fC")
		{
			bitMask = 0ll << shiftIndex;
		}
		else if ((*itChan)("Gain").getStringValue() == "240fC")
		{
			bitMask = 1ull << shiftIndex;
		}
		else if ((*itChan)("Gain").getStringValue() == "1pC")
		{
			bitMask = 2ull << shiftIndex;
		}
		else if ((*itChan)("Gain").getStringValue() == "10pC")
		{
			bitMask = 3ull << shiftIndex;
		}
	}
	return bitMask;
}

/**
 * analyse le fichier de parametres et recupere les valeurs correspndants
 * au registre 7 de aget
 */
uint64_t TestController::agetChanBitmask_reg7(int iChan,
		CCfg::CConfig::Iterator& itChan)
{
	uint64_t bitMask = 0ull;
	if (iChan >= 34 and iChan <= 67 and iChan != 45 and iChan != 56)
	{
		int shiftIndex = 2 * (iChan - 34) - 4;
		if (iChan < 45)
		{
			shiftIndex = 2 * (iChan - 34);
		}
		else if (iChan < 56)
		{
			shiftIndex = 2 * (iChan - 34) - 2;
		}

		if ((*itChan)("Gain").getStringValue() == "120fC")
		{
			bitMask = 0ll << shiftIndex;
		}
		else if ((*itChan)("Gain").getStringValue() == "240fC")
		{
			bitMask = 1ull << shiftIndex;
		}
		else if ((*itChan)("Gain").getStringValue() == "1pC")
		{
			bitMask = 2ull << shiftIndex;
		}
		else if ((*itChan)("Gain").getStringValue() == "10pC")
		{
			bitMask = 3ull << shiftIndex;
		}
	}
	return bitMask;
}
//__________________________________________________________________________________________________
/**
 * Builds mask for setting inhibition mode of channels 0 to 33 (except for FPN channels 11 and 22) in AGET register 10.
 * See section 9.7.K of AGET Data Sheet.
 * @param channelIndex Index of the channel to lookup in the configuration file.
 * @param channelCfg Channel configuration.
 * @return Mask where the 2 bits corresponding to the channel index have been set according to the configuration file.
 */
uint64_t TestController::agetChanBitmask_reg10(size_t channelIndex,
		CCfg::CConfig::Iterator& channelCfg)
{
	uint64_t reg10Value = UINT64_C(0);

	if (channelIndex > 32 or channelIndex == 11 or channelIndex == 22)
		return reg10Value;

	// Position for channel 2 bits in register
	size_t physicalChannelIndex = channelIndex;
	if (channelIndex > 11) physicalChannelIndex--;
	if (channelIndex > 22) physicalChannelIndex--;
	const size_t shiftIndex = 64-2*(physicalChannelIndex + 1);

	// Value of 2 bits
	// 0x0 = no inhibition
	// 0x1 = channel contributes neither to multiplicity nor to hit register
	// 0x2 = channel does not contribute to multiplicity
	// 0x3 is equivalent to 0x1
	std::string inhibitionMode = (*channelCfg)("TriggerInhibition").getStringValue();
	uint64_t modeValue = 0;
	if (inhibitionMode == "inhibit_channel" or inhibitionMode == "inhibit_trigger_function") // aliases for backward compatibility
	{
		modeValue = UINT64_C(0x1);
	}
	else if (inhibitionMode == "inhibit_trigger" or inhibitionMode == "inhibit_trigger_data") // aliases for backward compatibility
	{
		modeValue = UINT64_C(0x2);
	}

	// Set mask
	BitFieldHelper<uint64_t>::setField(reg10Value, shiftIndex, 2, modeValue);

	return reg10Value;
}
//__________________________________________________________________________________________________
/**
 * Builds mask for setting inhibition mode of channels 34 to 67 (except for FPN channels 45 and 56) in AGET register 11.
 * See section 9.7.L of AGET Data Sheet.
 * @param channelIndex Index of the channel to lookup in the configuration file.
 * @param channelCfg Channel configuration.
 * @return Mask where the 2 bits corresponding to the channel index have been set according to the configuration file.
 */
uint64_t TestController::agetChanBitmask_reg11(size_t channelIndex,
		CCfg::CConfig::Iterator& channelCfg)
{
	uint64_t reg11Value = UINT64_C(0);

	if (channelIndex < 34 or channelIndex > 67 or channelIndex == 45 or channelIndex == 56)
		return reg11Value;

	// Position for channel 2 bits in register
	size_t physicalChannelIndex = channelIndex;
	if (channelIndex > 45) physicalChannelIndex--;
	if (channelIndex > 56) physicalChannelIndex--;
	const size_t shiftIndex = 2*(physicalChannelIndex-34);

	// Value of 2 bits
	// 0x0 = no inhibition
	// 0x1 = channel contributes neither to multiplicity nor to hit register
	// 0x2 = channel does not contribute to multiplicity
	// 0x3 is equivalent to 0x1
	std::string inhibitionMode = (*channelCfg)("TriggerInhibition").getStringValue();
	uint64_t modeValue = 0;
	if (inhibitionMode == "inhibit_channel" or inhibitionMode == "inhibit_trigger_function") // aliases for backward compatibility
	{
		modeValue = UINT64_C(0x1);
	}
	else if (inhibitionMode == "inhibit_trigger" or inhibitionMode == "inhibit_trigger_data") // aliases for backward compatibility
	{
		modeValue = UINT64_C(0x2);
	}

	// Set mask
	BitFieldHelper<uint64_t>::setField(reg11Value, shiftIndex, 2, modeValue);

	return reg11Value;
}

/**
 * Extracts from configuration file the value corresponding to AGET register 12.
 */
uint64_t TestController::retrievalAget_reg12(int iCobo, int iAsad, int iAget)
{
	uint64_t datareg = 0ull;
	if (isAgetActive(iCobo, iAsad, iAget))
	{
		datareg = agetConfig(iCobo, iAsad, iAget)("Global")("ReadoutPointerOffset").getIntValue();
	}
	DBG_IF(verbose) << "Retrieved value '" << std::showbase << std::hex << datareg << std::dec << "' from config. for AGET register 'reg12'";
	return datareg;
}
//_________________________________________________________________________________________________
/**
 * Checks whether we should verify AGET slow control data have been correctly written.
 */
bool TestController::agetIsScCheck(int iCobo, int iAsad, int iAget)
{
	bool isCheck = agetConfig(iCobo, iAsad, iAget)("Control")("scCheckMode").getBoolValue();
	MSG_IF(verbose) << "Check of the configuration of AGET registers is " << (isCheck?"":"NOT ")<< "enabled";
	return isCheck;
}
//_________________________________________________________________________________________________
/**
 * Returns (initial) value of the difference between output voltage and reference voltage of the AsAd generator.
 *
 * @param coboIdx CoBo board identifier.
 * @return Value of the difference between output voltage and reference voltage of the AsAd generator [mV].
 */
float TestController::getAsadPulserVoltage(const size_t & coboIdx)
{
	float amplitude = 0;
	try
	{
		amplitude = coboConfig(coboIdx)("Generator")("amplitudeStart").getValueAsReal();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		LOG_WARN() << "Trying former configuration file format with AsAd no. 0";
		amplitude = asadConfig(coboIdx, 0)("Generator")("Amplitude", 0).getValueAsReal();
	}
	return amplitude;
}
//_________________________________________________________________________________________________
/**
 * Gets value of final amplitude for AsAd internal generator [mV].
 */
float TestController::getGeneAmpStop(const size_t & coboIdx)
{
	float amplitude = 0;
	try
	{
		amplitude = coboConfig(coboIdx)("Generator")("amplitudeStop").getValueAsReal();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		LOG_WARN() << "Trying former configuration file format with AsAd no. 0";
		amplitude = asadConfig(coboIdx, 0)("Generator")("AmplitudeStop").getValueAsReal();
	}
	return amplitude;
}
//_________________________________________________________________________________________________
/**
 * Gets value of amplitude increment for AsAd internal generator.
 * amplitude = amplitude + increment (increment and amplitude can be negative or positive).
 */
float TestController::getGeneIncr(const size_t & coboIdx)
{
	float amplitude = 0;
	try
	{
		amplitude = coboConfig(coboIdx)("Generator")("amplitudeStep").getValueAsReal();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		LOG_WARN() << "Trying former configuration file format with AsAd no. 0";
		amplitude = asadConfig(coboIdx, 0)("Generator")("Increment").getValueAsReal();
	}
	return amplitude;
}
//_________________________________________________________________________________________________
/**
 * Extracts from configuration the value to give to one of the readIfHit or readAlways CoBo registers.
 * 3 four-Byte CoBo registers in device 'ctrl': readIfHitMask1_4, readIfHitMask5_36 and readIfHit37_68
 * 3 four-Byte CoBo registers in device 'ctrl': readAlways1_4, readAlways5_36 and readAlways37_68
 *
 * @param filter Value of 'Reading' parameter to search for ('only_if_hit' or 'always').
 */
uint32_t TestController::getChannelReadingMask(const size_t & coboId, const size_t & asadId, const size_t & agetId,
		const size_t minChanId, const size_t & maxChanId, const std::string & filter)
{
	uint32_t regValue = 0;
	// Loop over configurations of concerned channels
	for (size_t chanId=minChanId; chanId <= maxChanId; ++chanId)
	{
		CCfg::CConfig cfg(agetConfig(coboId, asadId, agetId)("channel", chanId));
		if (cfg("isActive").getBoolValue())
		{
			const std::string chanParam = cfg("Reading").getStringValue();
			if (filter == chanParam)
			{
				BitFieldHelper<uint32_t>::setBit(regValue, maxChanId - chanId, true);
			}
		}
	}
	return regValue;
}
//__________________________________________________________________________________________________
/**
 * Sets ReadIfHit registers for a given AGET chip according to values in configuration file.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board the AGET belongs to (0 to 3).
 * @param agetIdx Index of the AGET chip within its AsAd board (0 to 3).
 */
void TestController::setReadIfHitMask(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx)
{
	const uint32_t value1_4 = getChannelReadingMask(coboIdx, asadIdx, agetIdx, 0, 3, "only_if_hit");
	const uint32_t value5_36 = getChannelReadingMask(coboIdx, asadIdx, agetIdx, 4, 35, "only_if_hit");
	const uint32_t value37_68 = getChannelReadingMask(coboIdx, asadIdx, agetIdx, 36, 67, "only_if_hit");

	eccClient().setReadIfHitMask(asadIdx, agetIdx, value1_4, value5_36, value37_68);
}
//__________________________________________________________________________________________________
/**
 * Sets ReadAlways registers for a given AGET chip according to values in configuration file.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board the AGET belongs to (0 to 3).
 * @param agetIdx Index of the AGET chip within its AsAd board (0 to 3).
 */
void TestController::setReadAlwaysMask(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx)
{
	const uint32_t value1_4 = getChannelReadingMask(coboIdx, asadIdx, agetIdx, 0, 3, "always");
	const uint32_t value5_36 = getChannelReadingMask(coboIdx, asadIdx, agetIdx, 4, 35, "always");
	const uint32_t value37_68 = getChannelReadingMask(coboIdx, asadIdx, agetIdx, 36, 67, "always");

	eccClient().setReadAlwaysMask(asadIdx, agetIdx, value1_4, value5_36, value37_68);
}
//__________________________________________________________________________________________________
/**
 * Sets thresholds of the 68 channels for zero suppression.
 * @param coboIdx CoBo ID
 * @param asadIdx AsAd ID
 * @param agetIdx AGET chip ID
 */
void TestController::setZeroSuppressionThresholds(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx)
{
	LOG_INFO() << "Setting zero suppression thresholds for AGET no. " << asadIdx*4+agetIdx;
	try
	{
		for (size_t channelId=0; channelId < 68; ++channelId)
		{
			const uint16_t threshold = agetConfig(coboIdx, asadIdx, agetIdx)("channel", channelId)("zeroSuppressionThreshold").getIntValue();
			eccClient().setCoBoZeroSuppressionThreshold(asadIdx, agetIdx, channelId, threshold);
		}
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what() << " Zero suppression thresholds have not been configured.";
	}
}
//__________________________________________________________________________________________________
/**
 * Writes and reads AGET register 1 through CoBo and AsAd.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param wchk Whether to check the register has been set to the expected value.
 */
bool TestController::setAget_reg1(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx, bool wchk)
{
	// Get AGET chip version to customize configuration between AGET and ASTRE
	uint16_t agetVersion = eccClient().agetVersion(asadIdx, agetIdx);
	// Get desired value
	const uint64_t value = retrievalAget_reg1(coboIdx, asadIdx, agetIdx, agetVersion);

	std::ostringstream oss;
	oss << "Setting register 'reg1' of AGET no. " << std::dec << 4*asadIdx+agetIdx  << " to " << std::showbase << std::hex << value;
	LOG_INFO() << oss.str();

	// Set register
	eccClient().agetRegWSC(asadIdx, agetIdx, "reg1", value);
	// Check register
	if (wchk)
	{
		const uint64_t resp = eccClient().agetRegRSC(asadIdx, agetIdx, "reg1");
		if (value != resp)
		{
			oss.str("");
			oss << "Check of slow control register 'reg1' failed for AGET no. " << std::dec << 4*asadIdx+agetIdx
					<< " (" << std::showbase << std::hex << resp << " instead of " << value << std::dec << ")\n";
			LOG_ERROR() << oss.str();
			return false;
		}
	}
	return true;
}
//__________________________________________________________________________________________________
/**
 * Writes and reads an AGET register through CoBo and AsAd.
 * @param registerName Name of the AGET register.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param wchk Whether to check the register has been set to the expected value.
 * @param value Value to set the register to.
 */
bool TestController::setAgetRegisterToValue(const std::string & registerName, const std::size_t & /* coboIdx */, const size_t & asadIdx, const size_t & agetIdx, bool wchk, const uint64_t & value)
{
	std::ostringstream oss;
	oss << "Setting register '" << registerName << "' of AGET no. " << std::dec << 4*asadIdx+agetIdx  << " to " << std::showbase << std::hex << value;
	LOG_INFO() << oss.str();

	eccClient().agetRegWSC(asadIdx, agetIdx, registerName, value);

	if (wchk)
	{
		uint64_t resp = eccClient().agetRegRSC(asadIdx, agetIdx, registerName);
		if (value != resp)
		{
			oss.str("");
			oss << "Check of slow control register '" << registerName << "' failed for AGET no. " << std::dec << 4*asadIdx+agetIdx
					<< " (" << std::showbase << std::hex << resp << " instead of " << value << std::dec << ")";
			// Handle special case of bit 21 of register 2 failure
			if ("reg2" == registerName)
			{
				uint64_t regXor =  value ^ resp;
				if (regXor == UINT64_C(0x200000))
				{
					LOG_WARN() << oss.str();
					return true;
				}
			}
			LOG_ERROR() << oss.str();
			return false;
		}
	}
	return true;
}
//__________________________________________________________________________________________________
/**
 * Sets an AGET register.
 * @param registerName Name of the AGET register.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param wchk Whether to check the register has been set to the expected value.
 */
bool TestController::setAgetRegister(const std::string & registerName, const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx, bool wchk)
{
	uint64_t value;
	if (registerName == "reg1")
	{
		// Get AGET chip version to customize configuration between AGET and ASTRE
		uint16_t agetVersion = eccClient().agetVersion(asadIdx, agetIdx);
		value = retrievalAget_reg1(coboIdx, asadIdx, agetIdx, agetVersion);
	}
	else if (registerName == "reg2")
		value = retrievalAget_reg2(coboIdx, asadIdx, agetIdx);
	else if (registerName == "reg12")
		value = retrievalAget_reg12(coboIdx, asadIdx, agetIdx);
	else
		value = retrievalChan_reg(registerName, coboIdx, asadIdx, agetIdx);

	return setAgetRegisterToValue(registerName, coboIdx, asadIdx, agetIdx, wchk, value);
}
//__________________________________________________________________________________________________
/**
 * Sets AGET register 12 and CoBo 'readDepth' register.
 * @param registerName Name of the AGET register.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param wchk Whether to check the register has been set to the expected value.
 */
bool TestController::setAget_reg12(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx,
		bool wchk)
{
	const uint16_t readoutPointerOffset = retrievalAget_reg12(coboIdx, asadIdx, agetIdx);

	// Set AGET register 12
	// This 16 bits register permits to reduce the number of SCA memory cells during the readout phase
	// by adding an offset (9 bits) on the readout pointer address.
	bool status = setAgetRegister("reg12", coboIdx, asadIdx, agetIdx, wchk);

	// Set CoBo readDepth register (10 bits): number of time buckets to read
	const uint16_t depth = is2pModeEnabled(coboIdx) ? 256 - readoutPointerOffset : 512 - readoutPointerOffset;
	eccClient().setCoBoReadDepth(depth);
	return status;
}
//__________________________________________________________________________________________________
/**
 * Writes and reads AGET register 2 through CoBo and AsAd.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param wchk Whether to check the register has been set to the expected value.
 */
bool TestController::setAget_reg2(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx, bool wchk)
{
	return setAgetRegister("reg2", coboIdx, asadIdx, agetIdx, wchk);
}
//__________________________________________________________________________________________________
/**
 * Writes and reads AGET registers 3 and 4 through CoBo and AsAd.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param wchk Whether to check the register has been set to the expected value.
 */
bool TestController::setAget_reg34(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx, bool wchk)
{
	bool success = setAgetRegister("reg3", coboIdx, asadIdx, agetIdx, wchk)
		and setAgetRegister("reg4", coboIdx, asadIdx, agetIdx, wchk);
	return success;
}
//__________________________________________________________________________________________________
/**
 * Writes and reads AGET registers 3 and 4 through CoBo and AsAd.
 * Verification occurs only if specified in the configuration file.
 */
bool TestController::setAgetDyn_reg3(int64_t val, int iAget, bool wchk)
{
	bool success = setAgetRegisterToValue("reg3", 0, 0, iAget, wchk, val & 0x3FFFFFFFFLLU)
			and setAgetRegisterToValue("reg4", 0, 0, iAget, wchk, 0LLU); // reset register 4 to zero
	return success;
}
//__________________________________________________________________________________________________
/**
 * Writes and reads AGET registers 3 and 4 through CoBo and AsAd.
 * Verification occurs only if specified in the configuration file.
 */
bool TestController::setAgetDyn_reg4(int64_t val, int iAget, bool wchk)
{
	bool success = setAgetRegisterToValue("reg3", 0, 0, iAget, wchk, 0LLU) // reset register 3 to zero
	 and setAgetRegisterToValue("reg4", 0, 0, iAget, wchk, val & 0x3FFFFFFFFLLU);
	return success;
}
//__________________________________________________________________________________________________
/**
 * Writes and reads AGET registers 6 and 7 through CoBo and AsAd.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param wchk Whether to check the register has been set to the expected value.
 */
bool TestController::setAget_reg67(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx, bool wchk)
{
	bool success = setAgetRegister("reg6", coboIdx, asadIdx, agetIdx, wchk)
	 and setAgetRegister("reg7", coboIdx, asadIdx, agetIdx, wchk);
	return success;
}
//__________________________________________________________________________________________________
/**
 * Writes and reads AGET registers 8 and 9 through CoBo and AsAd.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param wchk Whether to check the register has been set to the expected value.
 */
bool TestController::setAget_reg89(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx, bool wchk)
{
	bool success = setAgetRegister("reg8", coboIdx, asadIdx, agetIdx, wchk)
		and setAgetRegister("reg8msb", coboIdx, asadIdx, agetIdx, wchk)
		and setAgetRegister("reg9", coboIdx, asadIdx, agetIdx, wchk)
		and setAgetRegister("reg9msb", coboIdx, asadIdx, agetIdx, wchk);
	return success;
}
//__________________________________________________________________________________________________
/**
 * Writes and reads AGET registers 10 and 11 through CoBo and AsAd.
 * @param coboIdx Index of the CoBo board.
 * @param asadIdx Index of the AsAd board.
 * @param agetIdx Index of the AGET within its AsAd container.
 * @param wchk Whether to check the register has been set to the expected value.
 */
bool TestController::setAget_reg1011(const std::size_t & coboIdx, const std::size_t & asadIdx, const std::size_t & agetIdx, bool wchk)
{
	bool success = setAgetRegister("reg10", coboIdx, asadIdx, agetIdx, wchk)
			and setAgetRegister("reg11", coboIdx, asadIdx, agetIdx, wchk);
	return success;
}
//__________________________________________________________________________________________________
/**
 * Computes value to give to AsAd LCM2 register from configuraton parameters.
 */
uint64_t TestController::getAsAdLcm2Value(const size_t & coboIdx, const size_t & asadIdx)
{
	if (isCoboActive(coboIdx))
	{
		// Get CKR frequency
		float freqCKR_MHz = coboConfig(coboIdx)("Module")("readingClockFrequency").getRealValue(25.);

		if (isAsadActive(coboIdx, asadIdx))
		{
			CCfg::CConfig asadCfg(asadConfig(coboIdx, asadIdx));

			// Get ADC phase shift
			int32_t adcClockPhaseShift_deg = asadCfg("Clocking")("phaseShiftAdcClock").getIntValue(270);

			// Get ADC delay
			int32_t adcClockDelay_ps = asadCfg("Clocking")("delayAdcClock").getIntValue(5335);

			// Compute LCM2 value
			uint64_t lcm2Value = asad::utl::buildLcm2RegisterValue(freqCKR_MHz, adcClockPhaseShift_deg, adcClockDelay_ps);

			// Patch LCM2 value
			const bool lcm2DebugMode = asadCfg("Clocking")("Debug")("LCM2")("debugModeEnabled").getBoolValue(false);
			if (lcm2DebugMode)
			{
				::get::asad::utl::Lcm2Config lcm2Config;
				asad::utl::parseAsAdLcm2DebugConfig(asadCfg, lcm2Config);
				asad::utl::patchLcm2RegisterValue(lcm2Value, lcm2Config);
			}

			return lcm2Value;
		}
	}
	return UINT64_C(0);
}
//__________________________________________________________________________________________________
/**
 * Adjusts the value of the LCM1 register of all active AsAd boards using the configuration file.
 * @param freqCKW_MHz Frequency of the CKW sampling clock provided by CoBo.
 */
void TestController::asadSamplingClockInit(const size_t & coboIdx)
{
	// Get CKR frequency
	const float freqCKW_MHz = coboConfig(coboIdx)("Module")("writingClockFrequency").getRealValue();
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
	{
		if (isAsadActive(coboIdx, asadIdx))
		{
			// Compute LCM1 value
			float samplingClockDivider = 1; // WCKn = CKW
			try
			{
				samplingClockDivider = asadConfig(coboIdx, asadIdx)("Clocking")("samplingClockDivider").getRealValue();
			}
			catch (const CCfg::CConfig::NotFound & e)
			{
				LOG_WARN() << e.what();
			}
			uint64_t lcm1Value = asad::utl::buildLcm1RegisterValue(freqCKW_MHz, samplingClockDivider);

			// Patch LCM1 value
			CCfg::CConfig asadCfg(asadConfig(coboIdx, asadIdx));
			const bool lcm1DebugMode = asadCfg("Clocking")("Debug")("LCM1")("debugModeEnabled").getBoolValue(false);
			if (lcm1DebugMode)
			{
				::get::asad::utl::Lcm1Config lcm1Config;
				asad::utl::parseAsAdLcm1DebugConfig(asadCfg, lcm1Config);
				asad::utl::patchLcm1RegisterValue(lcm1Value, lcm1Config);
			}

			// Select AsAd
			selectAsAdForSlowControl(asadIdx);
			// Configure register LCM1 of AsAd Local Clock Manager device
			eccClient().setAsAdLcm1Register(lcm1Value);
		}
	}
	eccClient().setAsAdSlowControlMask(asadMask_);
}
//__________________________________________________________________________________________________
/**
 * Adjusts the value of the LCM2 register of all active AsAd boards using the configuration file.
 */
void TestController::asadReadingClockInit(const size_t & coboIdx)
{
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		if (isAsadActive(coboIdx, asadIdx))
		{
			// Get value
			uint64_t lcm2Value = getAsAdLcm2Value(coboIdx, asadIdx);
			// Select AsAd
			selectAsAdForSlowControl(asadIdx);
			// Set LCM2 register
			eccClient().setAsAdLcm2Register(lcm2Value);
		}
	}
}
//__________________________________________________________________________________________________
/**
 * Checks in config. whether CoBo should be in partial or full readout mode.
 */
bool TestController::isFullReadoutModeEnabled(const size_t & coboIdx)
{
	// Check readout mode of all chips
	size_t partialReadoutCount = 0;
	size_t activeChipCount = 0;
	for (size_t asadIdx=0; asadIdx < 4; ++asadIdx)
	{
		for (size_t agetIdx=0; agetIdx < 4; ++agetIdx)
		{
			if (not isAgetActive(coboIdx, asadIdx, agetIdx)) continue;

			activeChipCount++;
			if (not agetConfig(coboIdx, asadIdx, agetIdx)("Global")("Reg1")("isAllChannelRead").getBoolValue())
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
			msg << "All AGET chips of CoBo board '" << coboIdx
				<< "' should be in the same (full or partial) readout mode! Check the 'isAllChannelRead' parameters.";
			throw mdaq::utl::CmdException(msg.str());
		}
	}
	return fullReadoutModeEnabled;
}
//__________________________________________________________________________________________________
bool TestController::isFPNRead(int iCobo, int iAsad, int iAget)
{
	if (isAgetActive(iCobo, iAsad, iAget))
	{
		return agetConfig(iCobo, iAsad, iAget)("Global")("Reg1")("isFPNRead").getBoolValue();
	}
	return true;
}
//__________________________________________________________________________________________________
void TestController::enableFPNRead(int iCobo, int iAsad, int iAget)
{
	eccClient().enableFPNRead(isFPNRead(iCobo,iAsad,iAget));
}
//__________________________________________________________________________________________________
int32_t TestController::getEventNb(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("nbEvent").getIntValue();
	}
	return 1;
}
//_________________________________________________________________________________________________
/**
 * Returns the difference between default output voltage 1.102 mV.
 *
 * @param coboIdx CoBo board identifier.
 * @return Value of the difference between default output voltage and reference voltage of the AsAd generator [mV].
 */
float TestController::getAsadPulserDefaultVoltage(const size_t & coboIdx)
{
	float defaultVoltage = 0.0f;
	try
	{
		defaultVoltage = coboConfig(coboIdx)("Generator")("defaultVoltage").getValueAsReal();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	return defaultVoltage;
}
//__________________________________________________________________________________________________
/**
 * Checks mode to select for the AsAd generator in the TCM0 register.
 *  - Mode "test": a step voltage goes directly to the AGET test input.
 *  - Mode "calibration": a current pulse goes through a capacitor to the AGET calibration input.
 *
 * @see Section 1.4 of GET-AS-002-0005: "AsAd Test and Calibration Management".
 *
 * @return Returns true for calibration mode, false for test mode.
 * @param coboIdx CoBo board identifier
 * @param asadIdx AsAd board identifier
 */
bool TestController::getAsadPulserCalibMode(const size_t & coboIdx, const size_t & asadIdx)
{
	bool isCalibMode = false;
	try
	{
		isCalibMode = asadConfig(coboIdx, asadIdx)("Generator")("isCalibMode").getBoolValue();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		isCalibMode = coboConfig(coboIdx)("Generator")("isCalibMode").getBoolValue();
	}
	return isCalibMode;
}
//__________________________________________________________________________________________________
/**
 * Parses value to give to bit RNG (aka Rg) of register TCM0 of AsAd TCM device.
 * If T/C (aka Calib) = 1 (e.g. current pulse mode),
 *   RNG = 0 --> pulse created from  1pF AGET external capacitor.
 *   RNG = 1 --> pulse created from 11pF AGET external capacitor.
 *
 * The values proposed mistakenly used to be 10 nF (RNG=0) and 10 pF (RNG=1).
 *
 * @see Page 10 of "AsAd Slow Control & Registers Mapping", A. Rebii & J. Pibernat, May 2010.
 * @see Section 1.4 of GET-AS-002-0005: "AsAd Test and Calibration Management".
 * @see Table 8     of GET-AS-002-0010.
 * @see Section 8.1 of GET-QA-000-0005: "AGET Data Sheet".
 *
 * @return Returns true if RNG should be set (11pF capacitor), false otherwise (1pF).
 * @param coboIdx CoBo board identifier
 * @param asadIdx AsAd board identifier
 */
bool TestController::isAsAdPulserRng11pF(const size_t & coboIdx, const size_t & asadIdx)
{
	// Get user parameter value
	std::string capacitanceStr = "1pF";
	try
	{
		capacitanceStr = asadConfig(coboIdx, asadIdx)("Generator")("injectValue").getStringValue();

	}
	// For a time, parameter was set at the CoBo level, identical for all four AsAd boards
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		capacitanceStr = coboConfig(coboIdx)("Generator")("injectValue").getStringValue();
	}

	// Interpret string
	const bool use11pF = ("10pF" == capacitanceStr or "11pF" == capacitanceStr);
	return use11pF;
}
//__________________________________________________________________________________________________
/**
 * Checks if TCM device should be configured for step voltage requests from AsA external signals or from loading the TCM2 register.
 * @see GET-QA-000-0005
 * @see GET-QA-000-0003
 *
 * @return Returns true if TRG should be set (external trigger mode), false otherwise (slow control mode).
 * @param coboIdx CoBo board identifier
 * @param asadIdx AsAd board identifier
 */
bool TestController::getAsAdPulserTriggerMode(const size_t & coboIdx, const size_t & /* asadIdx */)
{
	// Get user parameter value
	bool externalMode = false;
	try
	{
		externalMode = coboConfig(coboIdx)("Generator")("externalTriggerMode").getBoolValue();
	}
	// Parameter added on March 28th, 2014
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	return externalMode;
}
//__________________________________________________________________________________________________
/**
 * Returns TCM device output step triggering delay (5 bits).
 * @see GET-QA-000-0005
 *
 * @return Delay expressed in number of WCK periods (0 to 31).
 * @param coboIdx CoBo board identifier
 * @param asadIdx AsAd board identifier
 */
uint8_t TestController::getAsAdPulserTriggerDelay(const size_t & coboIdx, const size_t & asadIdx)
{
	// Get user parameter value
	uint8_t delay_wckPeriods = 0;
	try
	{
		delay_wckPeriods = asadConfig(coboIdx, asadIdx)("Generator")("triggerDelay").getIntValue();
	}
	// Parameter added on March 28th, 2014
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	return delay_wckPeriods;
}
//__________________________________________________________________________________________________
/**
 * Selects mode of AsAd generator based on configuration file.
 *
 * @param coboIdx CoBo board identifier
 * @param asadIdx AsAd board identifier
 */
void TestController::asadPulserConfigure(const size_t & coboIdx, const size_t & asadIdx)
{
	// Check mode: either the pulser signal goes directly to the AGET test input (T) or it goes through a calibrated external capacitor (C).
	const bool TC = getAsadPulserCalibMode(coboIdx, asadIdx);

	// Check which of the 1 pF and 11 pF external capacitances to use.
	const bool RNG = isAsAdPulserRng11pF(coboIdx, asadIdx);

	// Check whether pulser should be triggered by slow control modification of TCM2 or by external signal from ISP1
	const bool TRG = getAsAdPulserTriggerMode(coboIdx, asadIdx);

	// Check triggering delay.
	const uint8_t DLY = getAsAdPulserTriggerDelay(coboIdx, asadIdx);

	// Check default voltage
	const float defaultVoltage_mV = getAsadPulserDefaultVoltage(coboIdx);

	// Check pulse voltage
	const float voltage_mV = getAsadPulserVoltage(coboIdx);

	eccClient().asadConfigureTcmDevice(TC, RNG, TRG, DLY, defaultVoltage_mV, voltage_mV);
}
//__________________________________________________________________________________________________
/**
 * Enables or disables the modification of the hit register before readout.
 */
void TestController::setWriteHitRegisterEnabled(int iCobo)
{
	const bool enabled = coboConfig(iCobo)("Module")("enableWriteHittedregister").getBoolValue();
	LOG_INFO() << "Modification of hit pattern before readout has been " << (enabled?"enabled":"disabled");
	eccClient().setWriteHitRegisterEnabled(enabled);
}
//__________________________________________________________________________________________________
/**
 * Checks in config. file the desired value for the multiplicity threshold.
 */
uint32_t TestController::getCoBoMultThres(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("multiplicityThreshold").getIntValue();
	}
	return 0;
}
//__________________________________________________________________________________________________
/**
 * Checks in config. file the desired value for the multiplicity threshold of the 2p mode second trigger.
 */
uint32_t TestController::getMultThres2p(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("multiplicityThreshold_2p").getIntValue();
	}
	return 0;
}
//__________________________________________________________________________________________________
/**
 * Checks in config. file the desired value for the width of the multiplicity trigger sliding window.
 */
uint32_t TestController::getCoBoSlidWin(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("multWindowSize").getIntValue();
	}
	return 10;
}
//__________________________________________________________________________________________________
/**
 * Sets multiplicity mode configuration parameters.
 */
void TestController::setMultiplicityMode(int iCobo)
{
	const uint32_t threshold = getCoBoMultThres(iCobo);
	MSG_IF(verbose) << "Setting multiplicity threshold to " << threshold;
	eccClient().setMultThreshold(threshold);
	const uint32_t windowWidth = getCoBoSlidWin(iCobo);
	MSG_IF(verbose) << "Setting width of multiplicity sliding window to " << windowWidth;
	eccClient().setSlidingWin(windowWidth);
}
//__________________________________________________________________________________________________
/**
 * Sets multiplicity 2p mode configuration parameters.
 */
void TestController::setMultThres2p(int iCobo)
{
	const uint32_t secondThreshold = getMultThres2p(iCobo);
	MSG_IF(verbose) << "Setting multiplicity threshold to " << secondThreshold;
	eccClient().setMultThreshold2p(secondThreshold);
}
//_________________________________________________________________________________________________
/**
 * Sets the CoBo in partial/full readout mode.
 * All AGET chips within a CoBo should be configured accordingly (see AGET register 1 bit 14).
 */
void TestController::setReadOutMode(const size_t & coboIdx)
{
	eccClient().setFullReadOutMode(isFullReadoutModeEnabled(coboIdx));
}
//_________________________________________________________________________________________________
/**
 * Sets the readout depth (number of time buckets).
 *
 * If the corresponding parameter is not found, the number of buckets to read is set to 512 or 256 instead.
 */
void TestController::setReadOutDepth(const size_t & coboIdx)
{
	const bool is2pEnabled = is2pModeEnabled(coboIdx);
	uint16_t readoutDepth = is2pEnabled ? 256 : 512 ;
	try
	{
		readoutDepth = coboConfig(coboIdx)("Module")("readoutDepth").getIntValue();
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what() << " Using a readout depth of " << readoutDepth << " bucket(s) instead.";
	}
	if (is2pEnabled and readoutDepth > 256)
	{
		LOG_ERROR() << readoutDepth << " is not a valid readout depth in 2p mode. Using 256 bucket(s) instead.";
		readoutDepth = 256;
	}
	else if (not is2pEnabled and readoutDepth > 512)
	{
		LOG_ERROR() << readoutDepth << " is not a valid readout depth. Using 512 bucket(s) instead.";
		readoutDepth = 512;
	}
	LOG_INFO() << "Setting readout depth to " << readoutDepth << " time bucket(s)";
	eccClient().setCoBoReadDepth(readoutDepth);
}
//_________________________________________________________________________________________________
/**
 * Sets the readout offset to indicate in the header of data frames.
 *
 * Sets the CoBo 'readDepth' register to the value given to the register 12 of the first active AGET chip.
 */
void TestController::setReadOutIndex(const size_t & coboIdx)
{
	size_t agetId = eccClient().firstActiveChip();
	const size_t asadId = agetId/4;
	agetId = agetId%4;
	uint16_t readoutOffset = agetConfig(coboIdx, asadId, agetId)("Global")("ReadoutPointerOffset").getIntValue();
	LOG_INFO() << "Setting readout offset to " << readoutOffset << " time bucket(s)";
	eccClient().setCoBoReadOffset(readoutOffset);
}
//_________________________________________________________________________________________________
/**
 * Enables or disables zero suppression and sets its polarity.
 *
 * If the corresponding parameter is not found, zero suppression is disabled.
 */
void TestController::setZeroSuppressionEnabled(const size_t & coboIdx)
{
	bool enableSuppression = false;
	bool zeroSuppressInvert = false;
	try
	{
		enableSuppression = coboConfig(coboIdx)("Module")("enableZeroSuppression").getBoolValue();
		zeroSuppressInvert  = coboConfig(coboIdx)("Module")("zeroSuppressionInverted").getBoolValue();
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	LOG_INFO() << (enableSuppression?"Enabling":"Disabling") << " zero suppression";
	if (enableSuppression and zeroSuppressInvert)
	{
		LOG_INFO() << "CoBo will suppress sample values ABOVE zero suppression thresholds.";
	}
	eccClient().setCoBoZeroSuppressionMode(enableSuppression, zeroSuppressInvert);
}
//_________________________________________________________________________________________________
/**
 * Switches to slow control mode.
 */
void TestController::switchToSlowControl(const size_t & /* coboIdx */)
{
	eccClient().switchToSlowControl();
}
//_________________________________________________________________________________________________
/**
 * Switches to fast control mode.
 */
void TestController::switchToFastControl(const size_t & /* coboIdx */)
{
	eccClient().switchToFastControl();
}
//_________________________________________________________________________________________________
/**
 * Injects charge using internal generator.
 */
void TestController::injectCharge(int64_t deltaV_mV, const size_t & /* coboIdx */)
{

	eccClient().injectCharge(deltaV_mV);
}
//_________________________________________________________________________________________________
/**
 * Extracts the trigger mode from the configuration file.
 */
std::string TestController::getTriggerMode(int iCobo)
{
	std::string mode = "periodically";
	if (isCoboActive(iCobo))
	{

		mode = coboConfig(iCobo)("Module")("triggerMode").getStringValue();
		MSG_IF(verbose) << "triggerMode: " << mode;
	}
	return mode;
}
//_________________________________________________________________________________________________
/**
 * Extracts the secondary trigger mode from the configuration file.
 */
std::string TestController::getSecondaryTriggerMode(int iCobo)
{
	std::string mode = "noTrigger";
	if (isCoboActive(iCobo))
	{
		try
		{
			try
			{
				mode = coboConfig(iCobo)("Module")("secondaryTriggerMode").getStringValue();
			}
			catch (const CCfg::CConfig::NotFound &)
			{
				mode = coboConfig(iCobo)("Module")("externalTriggerMode").getStringValue();
			}
		}
		catch (const CCfg::CConfig::NotFound &)
		{
			mode = getTriggerMode(iCobo);
			LOG_WARN() << "Using same mode for primary and secondary trigger logics.";
		}
		MSG_IF(verbose) << "secondaryTriggerMode: " << mode;
	}
	return mode;
}
//_________________________________________________________________________________________________
/*
 * Extracts the value of the ScwMultDelay register from the configuration file.
 */
int TestController::getScwMultDelay(const int coboIdx)
{
	uint64_t delay = 200;
	if (isCoboActive(coboIdx))
	{
		try
		{
			delay = coboConfig(coboIdx)("Module")("scwMultDelay").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound&)
		{ /* keep default value*/
		}
	}
	return delay;
}
//_________________________________________________________________________________________________
/*
 * Extracts the value of the ScwScrDelay register from the configuration file.
 */
int TestController::getScwScrDelay(const int coboIdx)
{
	uint64_t delay_10ns = 0;
	if (isCoboActive(coboIdx))
	{
		try
		{
			delay_10ns = coboConfig(coboIdx)("Module")("scwScrDelay").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound&)
		{
			try
			{
				delay_10ns = coboConfig(coboIdx)("Module")("writeReadDelay").getIntValue();
			}
			catch (const CCfg::CConfig::NotFound&)
			{
				// Keep default value
			}
		}
	}
	MSG_IF(verbose) << "ScwScrDelay: " << delay_10ns;
	return delay_10ns;
}
//_________________________________________________________________________________________________
/**
 * Sets CoBo trigger modes.
 */
void TestController::setTriggerModes(int iCobo)
{
	eccClient().setTriggerMode(getTriggerMode(iCobo));
	eccClient().setSecondaryTriggerMode(getSecondaryTriggerMode(iCobo));
}
//_________________________________________________________________________________________________
/**
 * Sets the delay from rising edge of SCW signal to start of multiplicity sliding window.
 */
void TestController::setScwMultDelay(const int coboIdx)
{
	eccClient().setScwMultDelay(getScwMultDelay(coboIdx));
}
//_________________________________________________________________________________________________
/**
 * Sets the delay between SCW low and SCR high.
 */
void TestController::setScwScrDelay(const int coboIdx)
{
	eccClient().setScwScrDelay(getScwScrDelay(coboIdx));
}
//_________________________________________________________________________________________________
/**
 * Extracts automatic trigger period from the configuration file.
 * @param coboIdx Index of the CoBo board.
 * @return Period in units of 10 ns.
 */
uint32_t TestController::getTriggerPeriod_10ns(const size_t & coboIdx)
{
	uint32_t period_10ns = 100000000u; // 1s
	if (isCoboActive(coboIdx))
	{
		try
		{
			period_10ns = coboConfig(coboIdx)("Module")("triggerPeriod_10ns").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound&)
		{
			period_10ns = 100000*coboConfig(coboIdx)("Module")("triggerPeriod").getIntValue(); // ms
		}
	}
	return period_10ns;
}
//_________________________________________________________________________________________________
/**
 * Sets the period of the trigger periodic mode for a CoBo target.
 * @param coboIdx Index of the CoBo board.
 */
void TestController::setTriggerPeriod(const size_t & coboIdx)
{
	eccClient().triggerPeriod_10ns(getTriggerPeriod_10ns(coboIdx));
}
//_________________________________________________________________________________________________
/**
 * Checks value of trigger delay.
 */
int64_t TestController::getTriggerDelay_10ns(const size_t & iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("triggerDelay").getIntValue();
	}
	return 800;
}
//_________________________________________________________________________________________________
/**
 * Sets value of delay between trigger event and actual freeze of the AGET memory.
 */
void TestController::setTriggerDelay_10ns(const size_t & coboIdx)
{
	eccClient().setTriggerDelay_10ns(getTriggerDelay_10ns(coboIdx));
}
//_________________________________________________________________________________________________
/**
 * Checks value of 2p mode second trigger delay.
 */
int64_t TestController::getTriggerDelay2p_10ns(const size_t & coboIdx)
{
	if (isCoboActive(coboIdx))
	{
		return coboConfig(coboIdx)("Module")("triggerDelay_2p").getIntValue();
	}
	return 800;
}
//_________________________________________________________________________________________________
/**
 * Sets value of delay between 2nd trigger event and actual freeze of the AGET memory.
 */
void TestController::setTriggerDelay2p_10ns(const size_t & coboIdx)
{
	eccClient().setTriggerDelay2p_10ns(getTriggerDelay2p_10ns(coboIdx));
}
//_________________________________________________________________________________________________
/**
 * recupere la valeur de retard du trigger
 */
int64_t TestController::getTriggerTimeOut2p_10ns(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("triggerTimeOut_2p").getIntValue();
	}
	return 800;
}
//_________________________________________________________________________________________________
/**
 * programme le retard du trigger
 */
void TestController::setTriggerTimeOut2p_10ns(int iCobo)
{
	eccClient().setTriggerTimeOut2p_10ns(getTriggerTimeOut2p_10ns(iCobo));
}
//_________________________________________________________________________________________________
/**
 * recupere la valeur de retard du trigger
 */
int64_t TestController::getTriggerDeadTime2p_10ns(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("triggerDeadTime_2p").getIntValue();
	}
	return 800;
}
//_________________________________________________________________________________________________
/**
 * programme le retard du trigger
 */
void TestController::setTriggerDeadTime2p_10ns(int iCobo)
{
	eccClient().setTriggerDeadTime2p_10ns(getTriggerDeadTime2p_10ns(iCobo));
}

/**
 * recupere la valeur de retard du trigger
 */
int64_t TestController::getTriggerTime2p_10ns(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("triggerTime_2p").getIntValue();
	}
	return 800;
}
//_________________________________________________________________________________________________
/**
 * programme le retard du trigger
 */
void TestController::setTriggerTime2p_10ns(int iCobo)
{
	eccClient().setTriggerTime2p_10ns(getTriggerTime2p_10ns(iCobo));
}
//_________________________________________________________________________________________________
/**
 * Checks whether to enable the 2p mode.
 */
bool TestController::is2pModeEnabled(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("enableMem2pMode").getBoolValue();
	}
	return false;
}
//_________________________________________________________________________________________________
/**
 * Checks in the configuration file whether to subtract the moving average from the AGET multiplicity signal.
 */
bool TestController::getCoBoMultSubtractAverage(const size_t & coboIdx)
{
	bool enabled = false;
	if (isCoboActive(coboIdx))
	{
		try
		{
			enabled = coboConfig(coboIdx)("Module")("multSubtractAverage").getBoolValue();
		}
		catch (const CCfg::CConfig::NotFound&)
		{
			// keep default value
		}
	}
	return enabled;
}

/**
 * Checks in the configuration file how many moving deviations to subtract from the AGET multiplicity signal.
 */
uint8_t TestController::getCoBoMultSubtractDevs(const size_t & coboIdx)
{
	uint8_t value = 0;
	if (isCoboActive(coboIdx))
	{
		try
		{
			value = coboConfig(coboIdx)("Module")("multSubtractDevs").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound&)
		{
			// keep default value
		}
	}
	return value;
}
//_________________________________________________________________________________________________
/**
 * Configures the subtraction of the moving average and deviation from the multiplicity signal
 *  in the CoBo multiplicity trigger mode.
 */
void TestController::setCoBoMultSubtractMode(const size_t & coboIdx)
{
	bool subtractAverage = getCoBoMultSubtractAverage(coboIdx);
	uint8_t subtractDevs = getCoBoMultSubtractDevs(coboIdx);
	eccClient().setCoBoMultSubtractMode(subtractAverage, subtractDevs);
}
//_________________________________________________________________________________________________
/**
 * Checks in the configuration file the value of the suppressMultiplicity field.
 * Looks at CoBo[*].AsAd[*].Aget[*].Global.suppressMultiplicity
 */
uint8_t TestController::getCoBoSuppressMultiplicity(const size_t & coboIdx)
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
			if (not isAsadActive(coboIdx, asadIdx)) continue;


			if (isAgetActive(coboIdx, asadIdx, agetIdx))
			{
				try
				{
					activeChipCount++;
					if (agetConfig(coboIdx, asadIdx, agetIdx)("Global")("suppressMultiplicity").getBoolValue())
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
				msg << "Cannot suppress multiplicity from AGET position " << agetIdx << " on one of CoBo '" << coboIdx
								<< "' 's AsAd board while keeping it on other AsAd boards. Check parameter 'suppressMultiplicity'";
				throw mdaq::utl::CmdException(msg.str());
			}
		}
		// Suppress multiplicity
		BitFieldHelper<uint8_t>::setField(suppressMultiplicity, agetIdx, 1, suppressChipPosition ? 1 : 0);
	}
	return suppressMultiplicity;
}
//_________________________________________________________________________________________________
/**
 * Configures the 4 bits of the CoBo 'pipeCtrl' register indicating which of the 4 chips to ignore in the multiplicity trigger mode.
 *
 * On Micro-TCA CoBo, the same setting applies to all AsAd boards.
 */
void TestController::setCoBoSuppressMultiplicity(const size_t & coboIdx)
{
	const uint8_t suppressMultiplicity = getCoBoSuppressMultiplicity(coboIdx);
	eccClient().setCoBoSuppressMultiplicity(suppressMultiplicity);
}
//_________________________________________________________________________________________________
/**
 * Sets 2p mode in CoBo.
 */
void TestController::enableMem2pMode(const size_t & coboIdx)
{
	eccClient().enableMem2pMode(is2pModeEnabled(coboIdx));
}
//_________________________________________________________________________________________________
/**
 * Configures the value of the 'dataSource' field in the header of the data frames.
 * This Data source ID typically refers to the detector or detector type from which the frame data have been extracted.
 * FIXME: We should ensure that all boards use the same dataSource value.
 */
void TestController::setCoBoDataSourceId(const size_t & coboId)
{
	if (isCoboActive(coboId))
	{
		uint8_t dataSourceId = 0;
		try
		{
			dataSourceId = coboConfig(coboId)("Module")("dataSource").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound & e)
		{
			LOG_WARN() << e.what();
		}
		eccClient().setDataSourceId(dataSourceId);
	}
}
//_________________________________________________________________________________________________
/**
 * Configures the value of the CoBo index in the header of the data frames.
 * FIXME: We should ensure that each board uses a different ID value.
 */
void TestController::setCoBoId(const size_t & coboId)
{
	if (isCoboActive(coboId))
	{
		uint8_t coboIdentifier = 0;
		try
		{
			coboIdentifier = coboConfig(coboId)("Module")("coboId").getIntValue();
		}
		catch (const CCfg::CConfig::NotFound & e)
		{
			LOG_WARN() << e.what();
		}
		eccClient().setCoBoId(coboIdentifier);
	}
}
//_________________________________________________________________________________________________
/**
 * Gets clock frequency for chip readout from configuration.
 */
float TestController::getCoBoReadingClockFrequency(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("readingClockFrequency").getRealValue();
	}
	return 25.0;
}
//__________________________________________________________________________________________________
/**
 * Sets CoBo clock frequency for AGET chip readout.
 */
void TestController::setCoBoReadingClockFrequency(int iCobo)
{
	eccClient().setCoBoReadingClockFrequency(
			getCoBoReadingClockFrequency(iCobo));
}
//__________________________________________________________________________________________________
/**
 * Checks value of delay for AEGT chip readout in configuration.
 */
int32_t TestController::getReadDataDelay(int iCobo)
{
	if (isCoboActive(iCobo))
	{
		return coboConfig(iCobo)("Module")("readDataDelay").getIntValue();
	}
	return 16;
}
//__________________________________________________________________________________________________
void TestController::saveAsicNb(int asicN)
{
	asicNum = asicN;
	LOG_DEBUG() << asicNum << " is the # of the component being tested";
}
//__________________________________________________________________________________________________
/**
 *  Writes manual read data delay into CoBo board.
 */
void TestController::setReadDataDelay(int iCobo)
{
	eccClient().setReadDataDelay(getReadDataDelay(iCobo));
}
//__________________________________________________________________________________________________
/**
 * Fetches sampling frequency in configuration file.
 */
float TestController::getWritingClockFrequency(const size_t & coboId, bool & pllConfigEnabled)
{
	if (isCoboActive(coboId))
	{
		try
		{
			pllConfigEnabled = not coboConfig(coboId)("Module")("skipPLLConfig").getBoolValue();
		}
		catch (const CCfg::CConfig::NotFound & e)
		{
			if (coboVersion_ >= 2u)
			{
				LOG_WARN() << e.what();
			}
			pllConfigEnabled = true;
		}
		return coboConfig(coboId)("Module")("writingClockFrequency").getRealValue();
	}
	return 100.0;
}
//__________________________________________________________________________________________________
/**
 * Sets CoBo writing clock frequency.
 */
void TestController::setCoBoWritingClockFrequency(const size_t & coboId)
{
	// Check whether to skip configuration of CoBo PLL device
	bool pllConfigEnabled = true;
	// Check sampling frequency
	float freq_MHz = getWritingClockFrequency(coboId, pllConfigEnabled);
	// Check advanced PLL parameters
	get::cobo::PllRegisterMap regs;
	if (coboVersion_ >= 2)
	{
		for (size_t regIndex=10; regIndex <= 13; ++regIndex)
		{
			regs[regIndex] = rc::CoBoConfigParser::getPllRegValue(coboConfig(coboId), regIndex);
		}
	}

	eccClient().setWritingClockFrequency(freq_MHz, pllConfigEnabled, regs);
}
//__________________________________________________________________________________________________
/**
 * Sets function of the 4 LEMO connectors on the CoBo front panel.
 */
void TestController::setCoBoLemoModes(const size_t & coboId)
{
	if (coboVersion_ < 2u) return;

	LOG_INFO() << "Configuring LEMO connectors";
	try
	{
		for (size_t lemoIndex=0; lemoIndex<4u; lemoIndex++)
		{
			uint8_t lemoMode = rc::CoBoConfigParser::getLemoConnectorMode(coboConfig(coboId), lemoIndex);
			eccClient().setCoBoLemoMode(lemoIndex, lemoMode);
		}
	}
	catch (const mdaq::utl::CmdException & e)
	{
		LOG_WARN() << e.what();
	}
}
//__________________________________________________________________________________________________
/**
 * Extracts from the configuration file the parameters of the circular buffer.
 */
uint32_t TestController::getCircularBufferParameter(const std::string & param,
		const size_t & coboId)
{
	if (isCoboActive(coboId))
	{
		return coboConfig(coboId)("CircularBuffer")(param).getHexValue();
	}
	return 0;
}
//__________________________________________________________________________________________________
/**
 * Writes into the CoBo board the parameters of the circular buffer.
 */
void TestController::setCircularBuffers(const size_t & coboIdx)
{
	// Set parameters shared by all 4 buffers
	uint32_t eventsPerInterrupt = 1u;
	try
	{
		eventsPerInterrupt = getCircularBufferParameter("interruptRate", coboIdx);
	}
	catch (const CCfg::CConfig::NotFound &) {}
	const uint32_t levelAlmostFullMemory = getCircularBufferParameter("levelAlmostFullMemory", coboIdx);
	const uint32_t levelAlmostEmptyMemory = getCircularBufferParameter("levelAlmostEmptyMemory", coboIdx);
	eccClient().setCircularBufferCommonParameters(eventsPerInterrupt, levelAlmostFullMemory, levelAlmostEmptyMemory);

	// Try to customize the size of each buffer
	try
	{
		for (size_t asadIdx=0; asadIdx <4u; ++asadIdx)
		{
			const uint32_t bufferStart = asadConfig(coboIdx, asadIdx)("CircularBuffer")("startOfMemory").getHexValue();
			const uint32_t bufferEnd = asadConfig(coboIdx, asadIdx)("CircularBuffer")("endOfMemory").getHexValue();
			eccClient().setAsAdCircularBuffer(asadIdx, bufferStart, bufferEnd);
		}
	}
	// Divide memory in equal parts
	catch (const CCfg::CConfig::NotFound &)
	{
		eccClient().setCircularBuffers(
			getCircularBufferParameter("startOfMemory", coboIdx),
			getCircularBufferParameter("endOfMemory", coboIdx));
	}
}
//_________________________________________________________________________________________________
/**
 * Parses the configuration to extract parameters for setting the AsAd inspection manager.
 */
bool TestController::getInspectionLinesParameters(std::string& moduleName,
		std::string & signalNameISP1, std::string & signalNameISP2, int iCobo,
		int iAsad)
{
	try
	{
		if (isAsadActive(iCobo, iAsad))
		{
			moduleName = asadConfig(iCobo, iAsad)("InspectionLines")("deviceToInspect").getStringValue();
			signalNameISP1 = asadConfig(iCobo, iAsad)("InspectionLines")(moduleName)("ISP1").getStringValue();
			signalNameISP2 = asadConfig(iCobo, iAsad)("InspectionLines")(moduleName)("ISP2").getStringValue();
			return true;
		}
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}

	moduleName = "ADC";
	signalNameISP1 = "SPI_CS";
	signalNameISP2 = "CAS";
	return false;
}
//_________________________________________________________________________________________________
/**
 * Configures the board(s).
 */
void TestController::configure()
{
	// Initialize state of run loop
	initStopLoop();
	// Create file for data if necessary
	if (not ControllerOptions::instance().isDataRouterExternal())
		createOutputFiles();

	// Reconnect to node in case another client of the ECC changed current node
	eccClient().connectNode(testSummary().targetEndpoint().toString());

	const std::size_t coboIdx = 0;
	// Loop over each of the 4 AsAd boards
	for (std::size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		// Skip inactive AsAd board
		if (not isAsadActive(coboIdx, asadIdx))
			continue;
		LOG_INFO() << "Configuring AsAd board no. " << asadIdx;

		// Select AsAd board for slow control
		selectAsAdForSlowControl(asadIdx);

		// Loop over AGET chips
		switchToSlowControl(coboIdx);
		for (std::size_t agetIdx = 0; agetIdx < 4u; ++agetIdx)
		{
			// Skip inactive AGET
			if (not isAgetActive(coboIdx, asadIdx, agetIdx))
				continue;
			LOG_INFO() << "Configuring registers of AGET no. " << (4*asadIdx+agetIdx);

			// Select AGET chip for slow control
			selectAgetForSlowControl(asadIdx, agetIdx);

			// Enable/disable readout of AGET FPN channels
			enableFPNRead(coboIdx, asadIdx, agetIdx);
			// Set AGET chip registers
			const bool scCheck = agetIsScCheck(coboIdx, asadIdx, agetIdx);
			bool agetRegConfigOk = true;
			agetRegConfigOk = agetRegConfigOk and setAget_reg1(coboIdx, asadIdx, agetIdx, scCheck);
			agetRegConfigOk = agetRegConfigOk and setAget_reg2(coboIdx, asadIdx, agetIdx, scCheck);
			agetRegConfigOk = agetRegConfigOk and setAget_reg34(coboIdx, asadIdx, agetIdx, scCheck);
			agetRegConfigOk = agetRegConfigOk and setAget_reg67(coboIdx, asadIdx, agetIdx, scCheck);
			agetRegConfigOk = agetRegConfigOk and setAget_reg89(coboIdx, asadIdx, agetIdx, scCheck);
			agetRegConfigOk = agetRegConfigOk and setAget_reg1011(coboIdx, asadIdx, agetIdx, scCheck);
			agetRegConfigOk = agetRegConfigOk and setAget_reg12(coboIdx, asadIdx, agetIdx, scCheck);
			if (not agetRegConfigOk)
			{
				std::ostringstream oss;
				oss << "Error configuring registers of AGET no. " << asadIdx*4+agetIdx;
				throw ::mdaq::utl::CmdException(oss.str());
			}
			// Set CoBo registers for AGET chip 'readAlways' and 'readIfHit' masks
			setReadIfHitMask(coboIdx, asadIdx, agetIdx);
			setReadAlwaysMask(coboIdx, asadIdx, agetIdx);
			// Set CoBo thresholds for zero suppression mode
			setZeroSuppressionThresholds(coboIdx, asadIdx, agetIdx);
		}
		switchToSlowControl(coboIdx);

		// Select AsAd board for slow control
		selectAsAdForSlowControl(asadIdx);

		// Set AsAd generator mode
		asadPulserConfigure(coboIdx, asadIdx);

		// Read AsAd voltages, currents, temperatures
		readAsAdVdd(asadIdx);
		readAsAdVdd1(asadIdx);
		readAsAdIdd(asadIdx);
		const int16_t Ti = readAsAdIntTemp(asadIdx);
		LOG_DEBUG() << "Read AsAd " << asadIdx << " internal temperature: " << Ti << ' ' << char(0xB0) << 'C';
		const int16_t Te = readAsAdExtTemp(asadIdx);
		LOG_DEBUG() << "Read AsAd " << asadIdx << " external temperature: " << Te << ' ' << char(0xB0) << 'C';

		// Set AsAd inspection manager
		setAsAdInspectionManager(coboIdx, asadIdx);
	}

	// Reset AsAd and AGET masks
	eccClient().setAsAdSlowControlMask(asadMask_);
	eccClient().setAgetSlowControlMask(agetMask_);

	// Enable/disable hit register modification
	setWriteHitRegisterEnabled(coboIdx);
	// Enable/disable CoBo test acquisition mode
	setAcqMode(coboIdx);
	// Set parameters of multiplicity trigger mode
	setMultiplicityMode(coboIdx);
	setMultThres2p(coboIdx);
	// Set multiplicity trigger moving average/deviation subtraction
	setCoBoMultSubtractMode(coboIdx);
	// Select which chips to ignore in the multiplicity trigger
	setCoBoSuppressMultiplicity(coboIdx);
	// Enable/disable 2p mode in CoBo
	enableMem2pMode(coboIdx);
	// Set data source ID
	setCoBoDataSourceId(coboIdx);
	// Set CoBo ID
	setCoBoId(coboIdx);
	// Set trigger period, delays, ...
	setTriggerPeriod(coboIdx);
	setTriggerDelay_10ns(coboIdx);
	setTriggerDeadTime2p_10ns(coboIdx);
	setTriggerDelay2p_10ns(coboIdx);
	setTriggerTime2p_10ns(coboIdx);
	setTriggerTimeOut2p_10ns(coboIdx);
	// Set circular buffer parameters
	setCircularBuffers(coboIdx);
	// Select trigger mode
	setTriggerModes(coboIdx);
	setScwMultDelay(coboIdx);
	setScwScrDelay(coboIdx);
	// Set CoBo readout mode
	setReadOutMode(coboIdx);
	// Set number of time buckets to readout
	setReadOutDepth(coboIdx);
	// Set readout offset in data frame header
	setReadOutIndex(coboIdx);
	// Enable zero suppression or not
	setZeroSuppressionEnabled(coboIdx);

	// Select which circular buffers to read
	eccClient().daqSetCircularBuffersEnabled();

	// Initialize 2p mode
	eccClient().aget2pInitNow();

	// Configure AsAd periodic pulser
	configureAsAdPeriodicPulser();

	// Make all CoBo boards send 'ready' status to Mutant modules.
	// Must be done before configuring Mutant modules as it affects Mutant B's MRDY status bit.
	eccClient().setCoBoInitDone(true);
}
//__________________________________________________________________________________________________
void TestController::readAsAdTemperatures()
{
	// Read temperature of all active AsAd boards
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		if (not isAsadActive(0, asadIdx)) continue;
		const int16_t Ti = readAsAdIntTemp(asadIdx);
		const int16_t Te = readAsAdExtTemp(asadIdx);
		LOG_INFO() << "Read temperatures of AsAd board no. " << asadIdx << ": Tint = " << Ti << ' ' << char(0xB0) << "C, Text = " << Te << ' ' << char(0xB0) << 'C';
	}
}
//__________________________________________________________________________________________________
/**
 * Main loop of the test controller thread.
 * It executes the various algorithms associated to the types of tests.
 */
void TestController::run()
{
	try
	{
		// Reconnect to node in case another client of the ECC changed current node
		eccClient().connectNode(testSummary().targetEndpoint().toString());

		// Read temperature of all active AsAd boards
		readAsAdTemperatures();

		const std::string testType = testSummary().type();
		LOG_INFO() << "Test type: '" << testType << '\'';

		emit runStarted();
		isDataMonitoringEnabled_ = ControllerOptions::instance().isDataDisplayEnabled();
		std::ostringstream I2S;
		QString message("");
		stopLoop = false;
		int32_t amp = 0;
		int32_t incr = 0;
		int32_t AmpStop = 0;

		// Parameters for test of peaking time
		uint16_t peakingTimeIndex = 0;

		uint64_t chSelect = 1;
		int32_t threshold = getInitGloablThresh(0, 0, 0);
		int32_t ThresResult[128][72];
		uint8_t indexThreshold;
		uint16_t nLoop;
		volatile bool isFrameReceived = false;
		int8_t nullByteCnt;
		nbEvent = 0;
		iInject = 0;
		bool ischSelectReg3 = true;
		switchToSlowControl();
		uint32_t statNb = getEventNb(0);
		bool isMinimalThreshold = false;
		if (testSummary().type() == "Threshold")
		{
			if (getTriggerMode(0) == "periodically")
			{
				isMinimalThreshold = true;
			}
		}
		LOG_INFO() << "isMinimalThreshold value : " << isMinimalThreshold;
		MSG_IF(verbose) << "Event statistics: " << statNb;

		// Ensure full frame is immediately sent if we are waiting for it
		eccClient().daqSetAlwaysFlushData("Default" != testType and "BaseLine" != testType);

		// Get initial generator amplitude
		amp = getAsadPulserVoltage(0);
		// Get increment for generator
		incr = getGeneIncr(0);
		// Get final amplitude
		AmpStop = getGeneAmpStop(0);
		MSG_IF(verbose) << "Generator amplitude: " << amp << " mV";

		if (testType == "PeakingTime")
		{
			setPeakingTimeIndex(peakingTimeIndex);
		}
		else if (testType == "AutoGain")
		{
			chSelect = 1;
			setAgetDyn_reg3(chSelect, 0, true);
		}
		else if (testType == "Threshold")
		{
			indexThreshold = 0;
			nullByteCnt = 0;
			nLoop = 0;
			for (uint8_t a = 0; a < 128; a++)
				for (uint8_t z = 0; z < 72; z++)
					ThresResult[a][z] = 0.0;
		}

		// Start
		if (testType != "Threshold")
		{
			daqStart(true);
		}

		// If required , start periodic AsAd pulser
		if (isAsAdPeriodicPulserEnabled())
		{
			// Start
			eccClient().ecc()->startAsAdPeriodicPulser();
		}

		startAcqWindow = true;
		while (not stopLoop)
		{
			setRunFinished(false);
			if (testType == "slowcontrol")
			{
			}
			else if (testType == "StdAcquisition")
			{
				while ((nbEvent < statNb or statNb <= 0) and not stopLoop)
				{
					injectCharge(amp);
					injectCharge(0);

					while (not stopLoop)
					{
						isFrameReceived = extractFrame();
						if (isFrameReceived) break;
					}
					if (isFrameReceived and isDataMonitoringEnabled_)
					{
						msleep(2000);
					}
				}
				if (statNb > 0 and nbEvent >= statNb)
				{
					nbEvent = 0;
					setStop();
				}
			}
			else if (testType == "Acq2pMode")
			{
				// Fake a cancelled 2p event
				/*injectCharge(amp);
				injectCharge(0);
				QThread::sleep(2);*/

				while ((nbEvent < 2*statNb or statNb <= 0) and not stopLoop)
				{
					// Injection for first SCA half-memory
					injectCharge(amp);
					injectCharge(0);

					// Injection for second SCA half-memory
					injectCharge(amp/5.);
					injectCharge(0);

					// Wait for first frame
					if (extractFrame())
					{
						if (isDataMonitoringEnabled_) msleep(2000);
					}
					// Wait for second frame
					if (extractFrame())
					{
						if (isDataMonitoringEnabled_) msleep(2000);
					}
				}
				if (statNb > 0 and nbEvent >= 2*statNb)
				{
					nbEvent = 0;
					setStop();
				}
			}
			else if (testType == "InjectTest")
			{
				for (iInject = 0; iInject < 100 and not stopLoop; ++iInject)
				{
					injectCharge(amp);
					injectCharge(0);

				}
				LOG_INFO() << "Injected event no. " << iInject;

				for (int32_t i = 0; i < 100 and not stopLoop; ++i)
				{
					while (not stopLoop)
					{
						isFrameReceived = extractFrame();
						if (isFrameReceived) break;
					}
				}
			}
			else if (testType == "BaseLine" or testType == "Default")
			{
				while (not stopLoop)
				{
					isFrameReceived = extractFrame();
					if (isFrameReceived) break;
				}
			}
			else if (testType == "Gain")
			{
				injectCharge(amp);
				injectCharge(0);
				if (isDataMonitoringEnabled_) msleep(2000);
				while (not stopLoop)
				{
					isFrameReceived = extractFrame();
					if (isFrameReceived) break;
				}
				if (nbEvent == statNb)
				{
					MSG_IF(verbose) << "Stopping to send data (" << nbEvent << " event(s))";
					nbEvent = 0;
					daqStart(false);
					switchToSlowControl();
					amp = amp + incr;
					MSG_IF(verbose) << "Amplitude: " << amp << " mV";
					switchToFastControl();
					daqStart(true);

				}
				if ((amp > AmpStop))
				{
					MSG_IF(verbose) << "Stopping to send data (" << amp << " mV)";
					setStop();
				}

			}
			else if (testType == "AutoGain")
			{
				injectCharge(amp);
				injectCharge(0);
				if (isDataMonitoringEnabled_) msleep(2000);
				while (not stopLoop)
				{
					isFrameReceived = extractFrame();
					if (isFrameReceived) break;
				}
				if (nbEvent == statNb)
				{
					MSG_IF(verbose) << "Stopping to send data (" << nbEvent << " event(s))";
					nbEvent = 0;
					daqStart(false);
					switchToSlowControl();
					amp = amp + incr;
					MSG_IF(verbose) << "Amplitude: " << amp << " mV";
					switchToFastControl();
					daqStart(true);

				}
				if ((amp > AmpStop))
				{
					daqStart(false);
					switchToSlowControl();
					amp = getAsadPulserVoltage(0);
					LOG_INFO() << "Starting amplitude: " << amp;
					chSelect = chSelect << 1;

					if (ischSelectReg3)
					{
						 setAgetDyn_reg3(chSelect,0,true);

						if ((chSelect > 8589934593LLU) and (chSelect
								< 17179869185LLU))
						{
							chSelect = 1;
							ischSelectReg3 = false;
						}
					}
					if (not ischSelectReg3)
					{
						setAgetDyn_reg4(chSelect,0,true);
					}

					switchToFastControl();
					daqStart(true);
					LOG_INFO() << "New channel: " << chSelect;
				}
				if (chSelect > 8589934592LLU)
				{
					chSelect = 1;
					setStop();
				}

			}
			else if (testType == "Threshold")
			{
				daqStart(false);
				switchToSlowControl();
				LOG_INFO() << "Threshold test: slow control mode";
				uint8_t thMSB = ((threshold & 112) >> 4);
				setGlobalRegThreshold(0, 0, thMSB);
				uint8_t thLSB = (threshold & 15);
				setChanRegThreshold(0, 0, thLSB);
				LOG_INFO() << "Test threshold value: " << std::hex << threshold << " (" << thMSB << " - " << thLSB << ')' << std::dec;
				switchToFastControl();
				daqStart(true);
				if (not isMinimalThreshold)
				{
					injectCharge(amp);
					injectCharge(0);
				}
				emit startTps();
				LOG_INFO() << "Send startTimer";
				if (isDataMonitoringEnabled_) msleep(2000);
				nLoop++;
				while (startAcqWindow)
				{
					if (stopLoop)
						break;

					isFrameReceived = extractFrame();
					if (isFrameReceived)
					{
						// TODO: To be confirmed by Frederic Druillole
						mfm::DynamicBitset hitPat_0 = getFrame->headerField("hitPat_0").value< mfm::DynamicBitset >();
						for (unsigned ch = 0; ch < 68; ch++) // instead of 72 channels in CoBo 2.1b format masks
						{
							if (hitPat_0[ch])
							{
								++ThresResult[indexThreshold][ch];
							}
						}
						if (not isMinimalThreshold)
						{
							injectCharge(amp);
							injectCharge(0);
						}
					}
					else
					{
						nullByteCnt++;
						LOG_INFO() << "nullByte Counter: " << nullByteCnt;
					}
				}
				LOG_INFO() << "startAcqWindow " << startAcqWindow << ": number of events = " << nbEvent;
				MSG_IF(verbose) << "Stopping to send data (" << nbEvent << " event(s))";

				daqStart(false);

				MSG_IF(verbose) << "Threshold value: " << threshold;

				if ((threshold > 126) or (nullByteCnt > 4))
				{
					MSG_IF(verbose) << "Stopping to send data (" << threshold << " threshold)";
					//writeThresholdTable(ThresResult);
					setStop();
					QString temp = "Threshold Result:\n";
					for (int32_t index = 0; index < 128; index++)
					{
						temp += QString("Threshold(%1):").arg(index);
						for (uint8_t ich = 0; ich < 72; ich++)
						{
							temp += QString("%1 ").arg(ThresResult[index][ich]);
						}
						temp += QString(tr("\n"));
					}
					emit logHd(temp);
				}
				else
				{
					if (nLoop >= statNb)
					{

						std::ostringstream oss;
						oss << "resp (" << (int) indexThreshold << "):";
						for (unsigned ch = 0; ch < 68; ch++)
						{
							oss << std::dec << ThresResult[indexThreshold][ch] << " ";
						}
						LOG_INFO() << oss.str();
						threshold += 1;
						indexThreshold++;
						nbEvent = 0;
						nLoop = 0;
						LOG_INFO() << "New test threshold value: " << std::hex << threshold;
					}
					startAcqWindow = true;
				}
			}
			else if (testType == "PeakingTime")
			{
				while (not stopLoop)
				{
					while (nbEvent < statNb and not stopLoop)
					{
						// Inject
						LOG_DEBUG() << "Injecting...";
						injectCharge(amp);
						injectCharge(0);

						// Wait for frame
						if (extractFrame())
						{
							// Display
							if (isDataMonitoringEnabled_) msleep(1000);
						}
					}
					LOG_DEBUG() << "Received " << nbEvent << " frame(s)";
					if (nbEvent >= statNb)
					{
						// Stop DAQ
						daqStart(false);
						nbEvent = 0;

						// Compute FWHM

						// Set peaking time
						switchToSlowControl();
						peakingTimeIndex++;
						if (peakingTimeIndex < 16)
						{
							setPeakingTimeIndex(peakingTimeIndex);
						}
						switchToFastControl();

						// Clear extra data
						dataBloc.clear();

						// Restart DAQ
						nbEvent = 0;
						daqStart(true);

						if (peakingTimeIndex >= 16)
						{
							// Stop DAQ
							daqStart(false);

							// Write to log file

							// Exit test
							setStop();
						}
					}
				}
			}
			else if (testType == "hittedChannel")
			{
				LOG_INFO() << "Hit amplitude: " << std::dec << amp;//recuperation de l'amplitude generateur
				injectCharge(amp);
				injectCharge(0);
				if (isDataMonitoringEnabled_) msleep(2000);
				while (not stopLoop)
				{
					isFrameReceived = extractFrame();
					if (isFrameReceived) break;
				}
				if (nbEvent == statNb)
				{
					if (isVerbose())
						std::cout << "Stop Sending Data: " << nbEvent
								<< std::endl;
					nbEvent = 0;
					daqStart(false);
					switchToSlowControl();
					//amp = amp + incr;
					LOG_INFO() << "Test amplitude: " << std::dec << amp;//recuperation de l'amplitude generateur
					chSelect = chSelect << 1;

					if (ischSelectReg3)
					{
						setAgetDyn_reg3(chSelect, 0, true);

						if ((chSelect > 8589934593LLU) and (chSelect
								< 17179869185LLU))
						{
							chSelect = 1;
							ischSelectReg3 = false;
						}
					}
					if (not ischSelectReg3)
					{
						setAgetDyn_reg4(chSelect, 0, true);
					}

					switchToFastControl();
					daqStart(true);
					LOG_INFO() << "New Channel: " << chSelect;
				}
				if (chSelect > 8589934592LLU)
				{
					chSelect = 1;
					setStop();
				}
			}

		}

		// Read temperature of all active AsAd boards
		readAsAdTemperatures();

		setRunFinished(true);
		dataBloc.clear();
		LOG_INFO() << "Leaving test controller thread.";
	}
	CATCH_ALL_EXCEPTIONS()
}
//__________________________________________________________________________________________________
/**
 * fonction message pour l'arret de l'acquisition:
 * on arrete le thread, on ferme les fichiers d'ecriture si ouvert
 * on stoppe l'envoi des donnees depuis la cible
 */
void TestController::setEndTimer()
{
	try
	{
		startAcqWindow = false;
		LOG_INFO() << "End of timer: " << nbEvent << " event(s)";
	}
	CATCH_ALL_EXCEPTIONS()
}
//__________________________________________________________________________________________________
/**
 * Stops acquisition.
 *   - stops controller thread
 *   - stops data emission from target
 */
void TestController::setStop()
{
	try
	{
		stopLoop = true;
		setRunFinished(true);
		startAcqWindow = false;
		daqStart(false);
		MSG_IF(verbose) << "stopLoop=" << stopLoop;
		emit runFinished();
	}
	CATCH_ALL_EXCEPTIONS()
}
//__________________________________________________________________________________________________
CCfg::CConfig & TestController::coboConfig(const size_t & coboId)
{
	return runCfg_(coboInstanceName_, coboId);
}
//__________________________________________________________________________________________________
CCfg::CConfig & TestController::asadConfig(const size_t & coboId, const size_t & asadId)
{
	return coboConfig(coboId)("AsAd", asadId);
}
//__________________________________________________________________________________________________
CCfg::CConfig & TestController::agetConfig(const size_t & coboId, const size_t & asadId, const size_t & agetId)
{
	return asadConfig(coboId, asadId)("Aget", agetId);
}
//__________________________________________________________________________________________________
bool TestController::isCoboActive(const size_t & coboId)
{
	return coboConfig(coboId)("isActive").getBoolValue();
}
//__________________________________________________________________________________________________
/**
 * Checks in the configuration whether to trigger the AsAd pulser periodically.
 */
bool TestController::isAsAdPeriodicPulserEnabled(const size_t & coboId)
{
	bool enabled = false;

	try
	{
		enabled = coboConfig(coboId)("Generator")("periodicModeEnabled").getBoolValue();
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}

	return enabled;
}
//__________________________________________________________________________________________________
/**
 * Configures AsAd periodic pulser.
 */
void TestController::configureAsAdPeriodicPulser(const size_t & coboId)
{
	if (not isAsAdPeriodicPulserEnabled(coboId)) return;

	// Get amplitude
	const int32_t amplitude_mV = getAsadPulserVoltage(coboId);
	// Get period
	size_t period_ms = 1000u;
	try
	{
		period_ms = coboConfig(coboId)("Generator")("period_ms").getIntValue();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	// Get amplitude of second pulse
	bool doublePulseEnabled = false;
	int32_t secondPulseAmplitude_mV = 0;
	try
	{
		doublePulseEnabled = coboConfig(coboId)("Generator")("doublePulseEnabled").getBoolValue();
		secondPulseAmplitude_mV = coboConfig(coboId)("Generator")("secondPulseAmplitude").getRealValue();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}

	// Set amplitude and period
	LOG_INFO() << "AsAd pulser will be triggered every " << period_ms << " ms with a voltage of Vref + " << amplitude_mV << " mV";
	eccClient().ecc()->configureAsAdPeriodicPulser(amplitude_mV, period_ms);

	// Set second pulse
	if (doublePulseEnabled)
	{
		LOG_INFO() << "AsAd pulser will send a second pulse with a voltage of Vref + " << secondPulseAmplitude_mV << " mV";
	}
	try
	{
		eccClient().ecc()->configureAsAdDoublePulserMode(doublePulseEnabled?secondPulseAmplitude_mV:0);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		if (doublePulseEnabled)
		{
			throw;
		}
	}

	// Check whether the pulser should be periodic or Poisson
	bool randomMode = false;
	try
	{
		randomMode = coboConfig(coboId)("Generator")("poissonModeEnabled").getBoolValue();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	try
	{
		eccClient().ecc()->setRandomAsAdPulserEnabled(randomMode);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		if (randomMode)
		{
			throw;
		}
	}
}
//__________________________________________________________________________________________________
bool TestController::isAsadActive(const size_t & coboId, const size_t & asadId)
{
	return isCoboActive(coboId) and asadConfig(coboId, asadId)("Control")("isActive").getBoolValue();
}
//__________________________________________________________________________________________________
bool TestController::isAgetActive(const size_t & coboId, const size_t & asadId, const size_t & agetId)
{
	return isAsadActive(coboId, asadId) and agetConfig(coboId, asadId, agetId)("Control")("isActive").getBoolValue();
}
//__________________________________________________________________________________________________
/**
 * Waits for new data and copies frame for processing.
 */
bool TestController::extractFrame()
{
	try
	{
		// Wait until new frame is available or until timeout
		DataBloc::RawFrame* rawFrame = dataBloc.waitNextFrame(1000);
		if (rawFrame == 0)
		{
			MSG_IF(verbose) << "Timeout";
			return false;
		}

		// Create local copy of frame
		typedef io::basic_array_source< mfm::Byte > ByteSource;
		ByteSource source(&(*rawFrame)[0], (size_t) rawFrame->size());
		io::stream< ByteSource > dataIn(source);
		getFrame = mfm::Frame::read(dataIn);

        // Clear frame from raw data stack
		dataBloc.popFrame(rawFrame);

		// Save raw frame
		if (isSaveRawFrameEnabled_)
		{
			writeFrameData(*getFrame.get());
			// Check file size
			const std::size_t fileSize_B = rawFrameFile.tellp();
			// Create new file if maximum file size has been reached
			if (fileSize_B > maxFileSize_MiB*0x100000)
				createRawFrameFile("", false);
		}

		// Skip topology frames
		if (getFrame->header().frameType() > 2u)
		{
			return false;
		}

		// Monitoring
		if (ControllerOptions::instance().isDataDisplayEnabled() and frameDisplay != 0)
		{
			frameDisplay->frameAvailable(*getFrame);
		}

		// Increment event counter
		nbEvent++;

		return true;
	}
	CATCH_ALL_EXCEPTIONS()
	return false;
}

#undef CATCH_ALL_EXCEPTIONS

