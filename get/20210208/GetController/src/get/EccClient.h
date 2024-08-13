/* =====================================================================================================================
 * EccClient.h
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

#include "get/GetEcc.h"
#include <utl/net/SocketAddress.h>
#include "get/cobo/LedManager.h"
using get::cobo::LedType;
using get::cobo::LedState;
#include <Ice/Ice.h>
#include <cstddef> // std::size_t

class EccClient
{
private:
	static bool verbose;
	get::GetEccPrx ecc_;
	utl::net::SocketAddress sockAddr;
	uint16_t agetMask_; ///< Mask of active AGET chips
	uint16_t asadMask_; ///< Mask of active AsAd boards
	uint8_t userCalibChipMask_; ///< Which AGET to use for serial delay calibration for each AsAd
	const std::string servantName;
	static const size_t NUM_KNOWN_AGET_VERSIONS = 3;
	static const uint16_t KNOWN_AGET_VERSIONS[NUM_KNOWN_AGET_VERSIONS];
public:
	EccClient(const std::string& ipDotDecimalAndPort, const std::string& servantName = "Ecc");
	virtual ~EccClient() {}
	get::GetEccPrx ecc();

	/* Constants */
	enum ModeADC {DDR, SDR};
	size_t firstActiveChip() const;
	void connectNode(const std::string & targetAddress);
	void loadHwDescription(const std::string& targetAddress, const std::string& hwPath);
	void daqConnect(const std::string & dataRouterAddress, const std::string & flowType);
	void daqDisconnect();
	static void setVerbose(bool verbosity) { verbose = verbosity; };
	static bool isVerbose() { return verbose; }
	void setActiveAGetMask(const uint16_t mask);
	void setAgetSlowControlMask(const uint16_t & mask);
	void selectAgetForRead(const size_t & asadId, const size_t & agetId);
	bool isAsAdBoardActive(const size_t & asadIdx) const;
	uint16_t agetVersion(size_t asadId, size_t agetId);
	void selectCalibChips(const uint8_t calibChipMask);
	void setActiveAsAdMask(const uint16_t mask);
	void setAsAdSlowControlMask(const uint16_t & mask);
	void selectAsAdForRead(const size_t & asadId);
	void selectAsAdForSlowControl(const size_t & asadIdx);
	uint32_t checkHardwareVersion();
	uint32_t checkFirmwareRelease();
	void pipeInit();
	void triggerInit();
	void circularBufferInit();
	void coboMultReadMode(bool enable);
	void asadPowerStatus(const size_t & asadIdx, bool & plg, bool & al, bool & pwu);
	uint32_t asadGetIdentifier();
	void asadSetBoardUserID(const uint8_t & userID);
	void asadInit(int time_ms=1000, const uint8_t & abortIfNoPowerSupply=0xF);
	void asadPowerOff();
	void agetInit();
	void calibrateAutoReadDataDelay();
	void asadSetAutoReadDataDelayEnabled(const size_t & asadId, bool enabled);
	void asadSetISPM(const std::size_t & coboIdx, const std::size_t & asadIdx,
			const std::string & moduleName, const std::string & isp1SignalName, const std::string & isp2SignalName);
	void asadAdcInit(ModeADC modeADC = DDR);
	/// @}

	/// @name CoBo
	/// @{
	void setCoBoInitDone(bool done);
	void setCoBoMutantMode(bool enabled);
	/// @}

	/// @name AsAd Monitoring
	/// @{
	void asadMonitoringReset();
	void asadDisableAutoProtection(const size_t & asadIdx);
	void asadSetAutoProtectionLevel(const uint8_t & interruptMask1, const uint8_t &);
	void asadSetMonitoringMode(const uint8_t & config1, const uint8_t & config2, const uint8_t & config3);
	void startMonitorAlarm();
	bool getAsAdAlarmFlag(const size_t & asadIdx);
	uint8_t getAsAdAlertStatus();
	void asadMonitoringStop();
	uint8_t readVdd(const std::size_t & asadIdx);
	uint8_t readIdd(const std::size_t & asadIdx);
	uint8_t readVdd1(const std::size_t & asadIdx);
	int8_t readIntTemp(const std::size_t & asadIdx);
	int8_t readExtTemp(const std::size_t & asadIdx);
	void wVddLimit(uint8_t low, uint8_t high);
	void wIddLimit(uint8_t low, uint8_t high);
	void wVdd1Limit(uint8_t low, uint8_t high);
	void wIntTempLimit(uint8_t low, uint8_t high);
	void wExtTempLimit(uint8_t low, uint8_t high);
	/// @}

	/// @name CoBo LEDs
	/// @{
	void clearAllLeds();
	void setLedState(const LedType & type, const size_t & asadId, const LedState & state);
	/// @}

	/// @name Calibration
	/// @{
	void asadCal();
	void coboSetCalibrationAgetChips(const uint16_t & agetMask);
	void agetCal();
	void printAgetChipVersionNumbers(bool check=true);
	void adcCal(ModeADC modeADC);
	/// @}

	/// @name AsAd
	/// @{
	void asadDacInit();
	void asadConfigureTcmDevice(const bool & TC, const bool & RNG, const bool & TRG, const uint8_t & DLY, const float & defaultVoltage_mV, const float & voltage_mV);
	/// @}

	/// @name configuration functions
	/// @{
	void acqTestMode(bool status);
	void triggerPeriod_10ns(const uint32_t & period_10ns);
	void setTriggerDelay_10ns(const uint32_t & triggerDelay_10ns);
	void setTriggerDelay2p_10ns(const uint32_t & triggerDelay);
	void setTriggerTimeOut2p_10ns(int triggerDelay);
	void setTriggerDeadTime2p_10ns(int triggerDelay);
	void setTriggerTime2p_10ns(int triggerDelay);
	void setCoBoMultSubtractMode(bool subtractAverage, const uint8_t & subtractDevs);
	void setCoBoSuppressMultiplicity(const uint8_t & suppressMultiplicity);
	void enableMem2pMode(bool valid2pMode);
	void setDataSourceId(const uint8_t & dataSourceId);
	void setCoBoId(const uint8_t & coboId);
	void enableFPNRead(bool isFPNRead);
	void setTriggerMode(const std::string & mode);
	void setSecondaryTriggerMode(const std::string & mode);
	void setScwMultDelay(const uint64_t & delay);
	void setScwScrDelay(const uint64_t & delay_10ns);
	void setFullReadOutMode(bool isAllChannelRead);
	void setCircularBuffers(const uint32_t & startOfMemory, const uint32_t & endOfMemory);
	void setAsAdCircularBuffer(const size_t & asadId, const uint32_t & startOfMemory, const uint32_t & endOfMemory);
	void setCircularBufferCommonParameters(const uint32_t & eventsPerInterrupt, const uint32_t & levelAlmostFullMemory,
			const uint32_t & levelAlmostEmptyMemory);

	/// @name functions
	/// @{
	void daqCtrl(bool status, bool resetTime = false);
	void daqSetAlwaysFlushData(bool enable=false);
	void daqSetCircularBuffersEnabled();
	void led(void);
	bool isAsAdConnected(const size_t & asadIdx);
	bool isAsAdPowerSupplied(const size_t & asadIdx);
	bool isAsAdPoweredOn(const size_t & asadIdx);
	void agetRegWSC(const std::size_t & asadIdx, const std::size_t & agetIdx, const std::string & registerName, uint64_t value);
	void agetRegWSC(const std::string & registerName, uint64_t value);
	uint64_t agetRegRSC(const std::size_t & asadIdx, const std::size_t & agetIdx, const std::string & registerName);
	void setWriteHitRegisterEnabled(bool enabled=true);
	void setReadIfHitMask(const std::size_t & asadIdx, const std::size_t & agetIdx, const uint32_t & val1_4, const uint32_t &  val5_36, const uint32_t &  val37_68);
	void setReadAlwaysMask(const std::size_t & asadIdx, const std::size_t & agetIdx, const uint32_t &  val1_4, const uint32_t &  val5_36, const uint32_t &  val37_68);
	void adcTest(ModeADC modeADC);
	void asadSetAgetInputProtections();
	void asadConfigureAgetInputs(const uint8_t & adc0, const uint8_t & IM0, const uint8_t & IM1);
	void adcReset();
	void setCoBoReadingClockFrequency(float ckr_MHz);
	void setAsAdLcm1Register(const uint64_t & value);
	void setAsAdLcm2Register(uint64_t value_LCM2);
	void setReadDataDelay(int delay);
	void setWritingClockFrequency(float ckw_MHz, bool pllConfigEnabled, const get::cobo::PllRegisterMap & regs);
	void setCoBoLemoMode(const size_t lemoIndex, const uint8_t & mode);
	void switchToFastControl();
	void switchToSlowControl();
	void injectCharge(long deltaV_mV);
	void triggerNow();
	void setMultThreshold(int thres);
	void setMultThreshold2p(int thres);
	void setSlidingWin(int slid);
	void setCoBoReadDepth(const uint16_t & depth);
	void setCoBoReadOffset(const uint16_t & offset);
	void setCoBoZeroSuppressionMode(const bool & enable, const bool & invert);
	void setCoBoZeroSuppressionThreshold(const size_t & asadIdx, const size_t & agetIdx, const size_t channelIdx, const uint16_t & threshold);
	void aget2pInitNow();
	/// @}

	/// @name miscellaneous
	/// @{
	void sleep_ms(int time_ms);
	/// @}
};
