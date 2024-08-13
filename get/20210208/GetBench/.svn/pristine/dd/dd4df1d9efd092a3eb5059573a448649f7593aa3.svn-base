/**
 * @file CoBoNode.h
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

#ifndef get_rc_CoBoNode_h_INCLUDED
#define get_rc_CoBoNode_h_INCLUDED

#include "get/rc/CoBoConfigParser.h"
#include "get/cobo/CtrlNode.h"
#include "get/cobo/MutantLinkManager.h"
#include "get/GetEcc.h"
#include "get/rc/NodeId.h"
#include "get/rc/Node.h"
#include "mdaq/hw/RemoteRegisterProcessor.h"
#include "utl/net/SocketAddress.h"
#include <boost/cstdint.hpp>
#include <boost/ptr_container/ptr_map.hpp>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
/**
 * A class to manage a CoBo board.
 */
class CoBoNode: public get::rc::Node, protected get::rc::CoBoConfigParser
{
public:
	static const std::string NODE_TYPE;
	typedef boost::ptr_map<std::string, CoBoNode> Map;

public:
	typedef enum {DDR, SDR} AdcMode;
	CoBoNode(const std::string& coboInstId, mdaq::rc::HardwareNode* _hwNode=0, get::rc::BackEnd* backend=0);
	virtual ~CoBoNode();

	DataSenderId dataSenderId() const;
	virtual DataSenderSet dataSenders() const;

	void runInit();

	virtual void daqConnect();
	void setCoBoWritingClockFrequency(const float & frequency_MHz, bool pllConfigEnabled=true, const cobo::PllRegisterMap & regs = cobo::PllRegisterMap());
	void setCoBoReadingClockFrequency(const float& frequency_MHz);
	void switchToFastControl(uint16_t agetMask);
	void switchToSlowControl(uint16_t agetMask);
	void switchToFastControl();
	void switchToSlowControl();
	void selectSingleAgetReadoutChannel(uint16_t agetMask, size_t chanId);
	void selectSingleAgetReadoutChannel(size_t chanId);
	void selectSingleAgetTestChannel(size_t chanId);
	void setAgetLsbThreshold(const size_t & chanId, const uint8_t threshold);
	void setInitDone(bool done);
	void setMutantMode(bool enabled);
	void setMutantLink(::get::cobo::MutantLinkMode mode);
	void setLevel2Trigger(bool enabled);
	friend class AmcSlotManager;

protected:
	virtual void impl_prepare(CCfg::CConfig nodeCfg);
	virtual void impl_unprepare();
	virtual void impl_configure(CCfg::CConfig nodeCfg);
	virtual void impl_start();
	virtual void impl_stop();

	void postStart();
	void preStop();

protected:
	/// @name CoBo
	/// @{
	void setCoBoCkwDivider(const float & frequency_MHz);
	void setLemoModes(CCfg::CConfig nodeCfg);
	void circularBufferInit();
	void setReadDataDelay(CCfg::CConfig nodeCfg);
	uint32_t readHardwareVersion();
	uint32_t checkFirmwareRelease();
	void checkReleaseNumbers();
	void pipeInit();
	void pipeInterruptEnable(bool enable);
	void triggerInit();
	void calibrateAutoReadDataDelay(CCfg::CConfig nodeCfg);
	void calibrateAutoReadDataDelay();
	void selectReadDataDelay(CCfg::CConfig nodeCfg);
	void asadSetAutoReadDataDelayEnabled(size_t asadId, bool);
	void coboEnableFPNRead(CCfg::CConfig & nodeCfg);
	void coboEnableFPNRead(bool isFPNRead);
	bool isFullReadoutModeEnabled(CCfg::CConfig & nodeCfg);
	void setReadoutMode(CCfg::CConfig & nodeCfg);
	void setReadoutDepth(CCfg::CConfig & nodeCfg);
	void setReadoutOffset(CCfg::CConfig & nodeCfg);
	void agetSetReadIfHitMask(CCfg::CConfig & nodeCfg, size_t asadId, size_t agetId);
	void agetSetReadAlwaysMask(CCfg::CConfig & nodeCfg, size_t asadId, size_t agetId);
	void setWriteHitRegisterMask(CCfg::CConfig nodeCfg);
	void setTestMode(CCfg::CConfig & nodeCfg);
	void setCoBoMultSubtractMode(CCfg::CConfig & nodeCfg);
	void setCoBoSuppressMultiplicity(CCfg::CConfig & nodeCfg);
	void setCoBoZeroSuppressionMode(CCfg::CConfig & nodeCfg);
	void setCoBoZeroSuppressionThresholds(CCfg::CConfig & nodeCfg, size_t asadIdx, size_t agetIdx);
	void enableMem2pMode(CCfg::CConfig & nodeCfg);
	void setCircularBuffers(CCfg::CConfig & nodeCfg);
	void setCircularBuffersEnabled(uint16_t mask);
	void resetCoBoReadPointers();
	void aget2pInitNow();
	void coboSetCalibrationAgetChips(uint16_t agetMask, CCfg::CConfig nodeCfg);
	void subscribeToCoBoAlarms(CCfg::CConfig & nodeCfg);
	void clearAllLeds();
	void modifyLed(const cobo::LedType & type, size_t asadId, const cobo::LedState & state);
	/// @}

	/// @name Data
	/// @{
	void setDataSourceId(CCfg::CConfig & nodeCfg);
	void setCoBoId(CCfg::CConfig & nodeCfg);
	/// @}

	/// @name Trigger
	/// @{
	void setTriggerPeriod(CCfg::CConfig & nodeCfg);
	void setTriggerDelay(CCfg::CConfig & nodeCfg);
	void setTriggerDelay2p(CCfg::CConfig & nodeCfg);
	void setTriggerDeadTime2p(CCfg::CConfig & nodeCfg);
	void setTriggerTime2p(CCfg::CConfig & nodeCfg);
	void setTriggerTimeOut2p(CCfg::CConfig & nodeCfg);
	void setTriggerModes(CCfg::CConfig & nodeCfg);
	void setTriggerMode(CCfg::CConfig & nodeCfg);
	void setSecondaryTriggerMode(CCfg::CConfig & nodeCfg);
	void setScwMultDelay(CCfg::CConfig & nodeCfg);
	void setScwScrDelay(CCfg::CConfig & nodeCfg);
	void setMultiplicityMode(CCfg::CConfig nodeCfg);
	void setMultiplicityThreshold2p(CCfg::CConfig nodeCfg);
	/// @}

	/// @name AsAd configuration
	/// @{
	void asadInit(CCfg::CConfig nodeCfg, size_t delay_ms=10);
	void asadCal();
	void asadDacInit();
	void clockInit(CCfg::CConfig nodeCfg);
	void asadSamplingClockInit(CCfg::CConfig nodeCfg, const float & freqCKW_Mhz);
	void asadReadingClockInit(CCfg::CConfig nodeCfg, float freqCKR_Mhz);
	void adcReset();
	void adcTest(AdcMode modeADC);
	void asadAdcInit(AdcMode modeADC);
	void adcCal(AdcMode adcMode);
	void asadSetAgetInputProtections();
	void asadAgetInputManagerInit(CCfg::CConfig nodeCfg);
	void asadInputsManagerConfig(const uint8_t & adc0, const uint8_t & IM0, const uint8_t & IM1);
	void asadPowerStatus(const size_t & asadIdx, bool & plg, bool & al, bool & pwu);
	uint32_t asadGetIdentifier();
	void asadSetBoardUserID(const uint8_t & userID);
	void configureAsAdTcmDevice(CCfg::CConfig & nodeCfg, size_t asadId);
	void setAsAdPeriodicPulser(CCfg::CConfig & nodeCfg);
	/// @}

	/// @name AsAd monitoring
	/// @{
	void asadMonitoringInit(CCfg::CConfig nodeCfg);
	void asadAlertsInit(CCfg::CConfig nodeCfg);
	void asadMonitoringRead();
	void asadMonitoringCheckVDD(CCfg::CConfig nodeCfg, const size_t maxReads = 10, const double delay_ms=100);
	void asadMonitoringReset();
	void asadMonitoringSetLimits(CCfg::CConfig & nodeCfg, size_t asadId);
	void setAsAdVddLimits(CCfg::CConfig nodeCfg, size_t asadId);
	void setAsAdVdd1Limits(CCfg::CConfig nodeCfg, size_t asadId);
	void setAsAdIddLimits(CCfg::CConfig nodeCfg, size_t asadId);
	void setAsAdIntTempLimits(CCfg::CConfig nodeCfg, size_t asadId);
	void setAsAdExtTempLimits(CCfg::CConfig nodeCfg, size_t asadId);
	void asadDisableAutoProtection(const size_t & asadId);
	void asadSetAutoProtectionLevel(CCfg::CConfig & nodeCfg, size_t asadId);
	void asadMonitoringSetOperationalMode(CCfg::CConfig & nodeCfg, size_t asadIdx);
	void asadMonitoringSetOperationalMode(bool enableAlarms=true, bool averagingDisabled=false, bool singleChannelMode=false, uint8_t singleChannelItem=0, bool fastMode=false);
	double readAsAdVdd(size_t asadId);
	double readAsAdVdd1(size_t asadId);
	double readAsAdIdd(size_t asadId);
	int8_t readAsAdIntTemp(size_t asadId);
	int8_t readAsAdExtTemp(size_t asadId);
	void readAsAdTemperatures();
	void setAsAdInspectionManager(CCfg::CConfig nodeCfg, const std::size_t & asadId);
	/// @}

	/// @name AGET
	/// @{
	void agetInit(CCfg::CConfig nodeCfg);
	void agetCal(CCfg::CConfig nodeCfg);
	void printAgetChipVersionNumbers(bool check=true);
	void setGlobalRegThresholds(const uint32_t & value);
	void setGlobalRegThreshold(size_t asadIdx, size_t agetIdx, const uint32_t & value);
	void agetSetRegister(size_t asadIdx, size_t agetIdx, const std::string & registerName, const uint64_t & value, bool check=true);
	void agetSetRegisters(CCfg::CConfig & nodeCfg, size_t asadId, size_t agetId, bool check=true);
	uint64_t agetRegRSC(size_t asadIdx, size_t agetIdx, const std::string & registerName);
	void agetRegWSC(size_t asadIdx, size_t agetIdx, const std::string & registerName, const uint64_t & value, const bool verbose=true);
	/// @}

	/// @name Management
	/// @{
	void setAgetMask(CCfg::CConfig nodeCfg);
	bool isAgetActive(size_t asadId, size_t agetId) const;
	uint16_t agetVersion(size_t asadId, size_t agetId);
	size_t firstActiveChip() const;
	bool isAsadActive(size_t asadId) const;
	bool isCoBoActive() const;
	bool isAsAdConnected(size_t asadIdx);
	uint8_t asadAlarmBit(const size_t & asadIdx);
	bool isAsAdPowerSupplied(size_t asadIdx);
	bool isAsAdPoweredOn(const size_t & asadIdx);
	void setAsAdSlowControlMask(uint16_t mask);
	void setAgetSlowControlMask(uint16_t mask);
	void selectAsAdForRead(size_t asadId);
	void selectAsAdForSlowControl(size_t asadIdx);
	void selectAgetForSlowControl(size_t asadIdx, const size_t agetIdx);
	/// @}
	struct ConfigData
	{
		bool asadPulserPeriodicModeEnabled;
		size_t asadPulserPeriod_ms;
		long asadPulserVoltage_mV;
		bool asadPulserDoubleModeEnabled;
		bool asadPulserRandomModeEnabled;
		long asadPulserVoltage2_mV;
		ConfigData();
	};

protected:
	get::cobo::CtrlNodePrx ctrlPrx();
	get::cobo::AsAdPulserMgrPrx asadPulserMgr();
	get::cobo::AsAdAlarmMonitorPrx asadAlarmMonitor();
	cobo::LedManagerPrx ledManager();
	cobo::MutantLinkManagerPrx mutantLinkManager();
	mdaq::hw::DevicePrx& ctrl();
	mdaq::hw::DevicePrx& asad();
	mdaq::hw::DevicePrx& aget();
	mdaq::hw::DevicePrx& zeroSuppressDevice(size_t asadIdx);
	static void sleep_ms(size_t delay_ms=1000);
	uint16_t agetMask_; ///< Mask of active AGET chips (16 bits)
	uint16_t asadMask_; ///< Mask of active AsAd boards (4 bits)
	bool isActive_;
	bool useExternalTimeReset_;
	ConfigData configData_;
	static const size_t NUM_KNOWN_AGET_VERSIONS = 3;
	static const uint16_t KNOWN_AGET_VERSIONS[NUM_KNOWN_AGET_VERSIONS];
};

//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
#endif /* get_rc_CoBoNode_h_INCLUDED */
