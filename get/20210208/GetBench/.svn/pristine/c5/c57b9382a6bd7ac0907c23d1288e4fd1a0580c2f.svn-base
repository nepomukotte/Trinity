/**
 * @file CoBoConfigParser.cpp
 * @date 16 oct. 2012
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

#include "get/rc/CoBoConfigParser.h"
#include "CCfg/Attribute.h"
#include "utl/Logging.h"
#include "utl/BitFieldHelper.hpp"
#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
namespace ba = boost::algorithm;
#include <iomanip>

using ::utl::BitFieldHelper;

inline uint32_t logtwo(uint32_t seed)
{
	uint32_t result = 0;
	int maxBits = 32;
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

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
CoBoConfigParser::CoBoConfigParser()
{
}
//_________________________________________________________________________________________________
bool CoBoConfigParser::verbose_ = false;
void CoBoConfigParser::setVerbose(bool enable)
{
	verbose_ = enable;
}
//_________________________________________________________________________________________________
CCfg::CConfig CoBoConfigParser::asadConfig(CCfg::CConfig cfg, const size_t & asadId)
{
	return cfg("AsAd", asadId);
}
//_________________________________________________________________________________________________
CCfg::CConfig CoBoConfigParser::agetConfig(CCfg::CConfig cfg, const size_t & asadId, const size_t & agetId)
{
	return asadConfig(cfg, asadId)("Aget", agetId);
}
//_________________________________________________________________________________________________
bool CoBoConfigParser::isAsadActive(const CCfg::CConfig& cfg, size_t asadId)
{
	return asadConfig(cfg, asadId)("Control")("isActive").getBoolValue();
}
//_________________________________________________________________________________________________
/**
 * Checks whether to use the automatically calibrated read data delays.
 *
 */
bool CoBoConfigParser::isSCRAutoDelayEnabled(CCfg::CConfig cfg, size_t asadId)
{
	bool enable = false;

	try
	{
		enable = asadConfig(cfg, asadId)("Control")("enableAutoReadDataDelay").getBoolValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		try
		{
			enable = cfg("Module")("enableAutoReadDataDelay").getBoolValue();
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
 * Convert channel ID (from 0 to 67 including FPN channels) to physical channel ID (from 0 to 63).
 * @param channelId Channel ID
 * @return Physical channel ID or -1 for FPN channel.
 */
size_t CoBoConfigParser::agetPhysicalChannelId(const size_t channelId)
{
	size_t physicalId = channelId;
	if (channelId == 11 or channelId == 22 or channelId == 45 or channelId == 56)
		physicalId = -1;
	else if (channelId > 56)
		physicalId -= 4;
	else if (channelId > 45)
		physicalId -= 3;
	else if (channelId > 22)
		physicalId -= 2;
	else if (channelId > 11)
		physicalId -= 1;
	return physicalId;
}
//_________________________________________________________________________________________________
/** Return AGET or ASTRE register 1 state1<3:6> field for given peaking time. */
uint16_t CoBoConfigParser::peakingTimeIndex(size_t peakingTime, uint16_t agetVersion)
{
	// "Time". Sets the peaking time (ns) of the shaper by switching resistors on the PZC  SK filters.
	MSG_IF(verbose_) << "Shaper peaking time: " << peakingTime << " ns";

	uint16_t peakingTimeIndex = 0;

	// ASTRE
	if (agetVersion == 0x203)
	{
		const size_t astre_times[16] = {68, 116, 180, 225, 343, 442, 541, 641, 739, 836, 943, 1039, 2034, 3985, 5770, 7700};
		const size_t* found = std::find(astre_times, astre_times+16, peakingTime);
		if (found == astre_times+16)
		{
			peakingTimeIndex = 6;
			LOG_WARN() << peakingTime << " ns  is not a valid ASTRE peaking time value! Using " << astre_times[peakingTimeIndex] << " ns instead.";
		}
		else
		{
			peakingTimeIndex = std::distance(astre_times, found);
		}
	}
	// AGET
	else
	{
		const size_t aget_times[16] = {70, 117, 232, 280, 334, 383, 502, 541, 568, 632, 721, 760, 831, 870, 976, 1014};
		const size_t* found = std::find(aget_times, aget_times+16, peakingTime);
		if (found == aget_times+16)
		{
			peakingTimeIndex = 8;
			LOG_WARN() << peakingTime << " ns  is not a valid AGET peaking time value! Using " << aget_times[peakingTimeIndex] << " ns instead.";

		}
		else
		{
			peakingTimeIndex = std::distance(aget_times, found);
		}
	}
	return peakingTimeIndex;
}
//_________________________________________________________________________________________________
/**
 * Parses a configuration file and computes the desired value for AGET 32-bit register 1.
 */
uint32_t CoBoConfigParser::agetReg1Value(CCfg::CConfig & cfg, const size_t & asadId, const size_t & agetId, const uint16_t & agetVersion)
{
	uint32_t regValue = UINT32_C(0);
	uint32_t fieldValue;
	CCfg::CConfig agetCfg = agetConfig(cfg, asadId, agetId);

	// "Icsa". If set, the nominal CSA bias current is multiplied by two.
	const bool isIcsaSet = agetCfg("Global")("Reg1")("isIcsax2").getBoolValue();
	MSG_IF(verbose_) << "Nominal CSA bias current doubled? " << isIcsaSet;
	fieldValue = isIcsaSet?0x1:0x0;
	BitFieldHelper<uint32_t>::setField(regValue, 0, 1, fieldValue);

	 // "Gain". LSb and MSb of the internal test capacitor in Test mode.
	const std::string testCapacitorValue = agetCfg("Global")("Reg1")("TestModeRange").getStringValue();
	// We expect Farad sub-units or, in old config files, Coulomb sub-units
	MSG_IF(verbose_) << "Internal test capacitor ('TestModeRange'): " << testCapacitorValue;
	if (ba::starts_with(testCapacitorValue, "120f") or ba::starts_with(testCapacitorValue, "120 f"))
		fieldValue = UINT32_C(0x0);
	else if (ba::starts_with(testCapacitorValue, "240f") or ba::starts_with(testCapacitorValue, "240 f"))
		fieldValue = UINT32_C(0x1);
	else if (ba::starts_with(testCapacitorValue, "1p") or ba::starts_with(testCapacitorValue, "1 p") )
		fieldValue = UINT32_C(0x2);
	else if (ba::starts_with(testCapacitorValue, "10p") or ba::starts_with(testCapacitorValue, "10 p"))
		fieldValue = UINT32_C(0x3);
	else
	{
		LOG_WARN() << "Unexpected value for internal test capacitor: " << testCapacitorValue << ". Using 120 fF instead.";
		fieldValue = UINT32_C(0x0);
	}
	BitFieldHelper<uint32_t>::setField(regValue, 1, 2, fieldValue);

	// "Time". Sets the peaking time (ns) of the shaper by switching resistors on the PZC  SK filters.
	const int peakingTime = agetCfg("Global")("Reg1")("peackingTime").getIntValue();
	fieldValue = peakingTimeIndex(peakingTime, agetVersion);
	BitFieldHelper<uint32_t>::setField(regValue, 3, 4, fieldValue);

	// "Test". Defines the test mode: none (0), calibration (01), test (10), functionality (11).
	const std::string testMode = agetCfg("Global")("Reg1")("TestModeSelection").getStringValue();
	MSG_IF(verbose_) << "Test mode: " << testMode;
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
	BitFieldHelper<uint32_t>::setField(regValue, 7, 2, fieldValue);

	// "Integrator_mode". Bit integrator mode
	const std::string integratorMode = agetCfg("Global")("Reg1")("IntegrationMode").getStringValue();
	MSG_IF(verbose_) << "IntegrationMode: " << integratorMode;
	if (integratorMode == "variable")
		fieldValue = 0x0;
	else
		fieldValue = 0x1;
	BitFieldHelper<uint32_t>::setField(regValue, 9, 1, fieldValue);

	// "SCApointer". Bit SCA pointer 0 to pointer 1.
	const std::string scaPointer = agetCfg("Global")("Reg1")("SCA_Pointer").getStringValue();
	MSG_IF(verbose_) << "SCA_Pointer: " << scaPointer;
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
	BitFieldHelper<uint32_t>::setField(regValue, 10, 2, fieldValue);

	// "SCAsplitting". Bit SCA splitting; if set, enables the 2p mode by dividing the memory in two blocks.
	if (agetCfg("Global")("Reg1")("SCA_Splitting").getBoolValue())
		fieldValue = 0x1;
	else
		fieldValue = 0x0;
	BitFieldHelper<uint32_t>::setField(regValue, 12, 1, fieldValue);
	MSG_IF(verbose_) << "SCA_Splitting: " << fieldValue;


	// "Mode32Channels". Bit 32 channels mode. If set, only the first 32 channels are used.
	if (agetCfg("Global")("Reg1")("is32channels").getBoolValue())
		fieldValue = 0x1;
	else
		fieldValue = 0x0;
	BitFieldHelper<uint32_t>::setField(regValue, 13, 1, fieldValue);
	MSG_IF(verbose_) << "is32channels";

	// "Readout_mode". Bit readout mode; switches between full and partial readout modes. If set, full readout.
	if (agetCfg("Global")("Reg1")("isAllChannelRead").getBoolValue())
		fieldValue = 0x1;
	else
		fieldValue = 0x0;
	BitFieldHelper<uint32_t>::setField(regValue, 14, 1, fieldValue);
	MSG_IF(verbose_) << "isAllChannelRead";

	// "FPNreadout". Bit FPN readout: if set, forces the readout of the FPN channels (e.g. even in partial readout mode).
	if (agetCfg("Global")("Reg1")("isFPNRead").getBoolValue())
		fieldValue = 0x1;
	else
		fieldValue = 0x0;
	BitFieldHelper<uint32_t>::setField(regValue, 15, 1, fieldValue);
	MSG_IF(verbose_) << "isFPNRead: " << fieldValue;

	// "Polarity". Bit polarity. Controls the value of the DC voltage levels. Usually set to zero?
	if (agetCfg("Global")("Reg1")("isPositivePolarity").getBoolValue())
		fieldValue = 0x1;
	else
		fieldValue = 0x0;
	BitFieldHelper<uint32_t>::setField(regValue, 16, 1, fieldValue);
	MSG_IF(verbose_) << "isPositivePolarity: " << fieldValue;

	// "Vicm". Input common mode voltage (VICM) of the analog output buffer. Usually 0 e.g. 1.25 V?
	const std::string vicmValue = agetCfg("Global")("Reg1")("Vicm").getStringValue();
	MSG_IF(verbose_) << "Vicm: " << vicmValue;
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
	BitFieldHelper<uint32_t>::setField(regValue, 17, 2, fieldValue);

	// "DACthreshold".
	MSG_IF(verbose_) << "GlobalThresholdValue";
	BitFieldHelper<uint32_t>::setField(regValue, 19, 3, agetCfg("Global")("Reg1")(
			"GlobalThresholdValue").getIntValue());

	// "DACpolarity".
	MSG_IF(verbose_) << "isThresholdSignedPositive";
	fieldValue = agetCfg("Global")("Reg1")("isThresholdSignedPositive")
			.getBoolValue()?0x1:0x0;
	BitFieldHelper<uint32_t>::setField(regValue, 22, 1, fieldValue);

	// "Trigger_veto". Trigger bit: specifies the use of the veto on the trigger building. Defines how long after a hit the channel can be read again.
	// Values: 0 (none), 01 ( 4 microseconds), 10 (hit register width).
	const std::string triggerVetoOption = agetCfg("Global")("Reg1")(
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
	MSG_IF(verbose_) << "TriggerVetoOption: " << triggerVetoOption;
	BitFieldHelper<uint32_t>::setField(regValue, 23, 2, fieldValue);

	// "Synchro_discri".
	fieldValue = agetCfg("Global")("Reg1")("isckwriteSynchro")
			.getBoolValue()?0x1:0x0;
	MSG_IF(verbose_) << "isCKWriteSynchro: " << fieldValue;
	BitFieldHelper<uint32_t>::setField(regValue, 25, 1, fieldValue);

	// "tot".
	fieldValue = agetCfg("Global")("Reg1")("isTOTActive").getBoolValue()? 0x1:0x0;
	BitFieldHelper<uint32_t>::setField(regValue, 26, 1, fieldValue);
	MSG_IF(verbose_) << "isTOTActive: " << fieldValue;

	// "Range_trigg_width"
	fieldValue = agetCfg("Global")("Reg1")("isTriggerWidth200ns").getBoolValue()?0x1:0x0;
	BitFieldHelper<uint32_t>::setField(regValue, 27, 1, fieldValue);
	MSG_IF(verbose_) << "isTriggerWidth200ns: " << fieldValue;

	// "trigg_width"
	const std::string triggerWidthRange = agetCfg("Global")("Reg1")("TriggerWidthRange")
			.getStringValue();
	MSG_IF(verbose_) << "TriggerWidthRange: " << triggerWidthRange;
	if (triggerWidthRange == "width1")
	{
		// "Lsb_trigg_width"
		BitFieldHelper<uint32_t>::setField(regValue, 28, 1, 0x1);
		// "Msb_trigg_width"
		BitFieldHelper<uint32_t>::setField(regValue, 29, 1, 0x0);
	}
	else if (triggerWidthRange == "width2")
	{
		// "Lsb_trigg_width"
		BitFieldHelper<uint32_t>::setField(regValue, 28, 1, 0x0);
		// "Msb_trigg_width"
		BitFieldHelper<uint32_t>::setField(regValue, 29, 1, 0x1);
	}
	else if (triggerWidthRange == "width3")
	{
		// "Lsb_trigg_width"
		BitFieldHelper<uint32_t>::setField(regValue, 28, 1, 0x1);
		// "Msb_trigg_width"
		BitFieldHelper<uint32_t>::setField(regValue, 29, 1, 0x1);
	}
	else // "width0"
	{
		// "Lsb_trigg_width"
		BitFieldHelper<uint32_t>::setField(regValue, 28, 1, 0x0);
		// "Msb_trigg_width"
		BitFieldHelper<uint32_t>::setField(regValue, 29, 1, 0x0);
	}

	// "External".
	const std::string externalLink = agetCfg("Global")("Reg1")("ExternalLink").getStringValue();
	MSG_IF(verbose_) << "ExternalLink: " << externalLink;
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
	BitFieldHelper<uint32_t>::setField(regValue, 30, 2, fieldValue);

	return regValue;
}
//__________________________________________________________________________________________________
/**
 * Parses a configuration file and computes the desired value for AGET 32-bit register 2.
 */
uint32_t CoBoConfigParser::agetReg2Value(CCfg::CConfig & cfg, const size_t & asadId,
		const size_t & agetId)
{
	uint32_t regValue = 0;
	uint32_t fieldValue = 0;
	CCfg::CConfig agetCfg = agetConfig(cfg, asadId, agetId);

	// "debug0" to debug2". These 3 bits select the internal signal to be view on scoop.
	const std::string dbgSelection = agetCfg("Global")("Reg2")("DebugSelection").getStringValue();
	MSG_IF(verbose_) << "DebugSelection:" << dbgSelection;
	if (dbgSelection == "standby")
		fieldValue = 0x0;
	else if (dbgSelection == "CSA")
		fieldValue = 0x1;
	else if (dbgSelection == "CR")
		fieldValue = 0x2;
	else if (dbgSelection == "Gain-2")
		fieldValue = 0x3;
	else if (dbgSelection == "PositiveInputDiscri")
		fieldValue = 0x4;
	else if (dbgSelection == "NegativeInputDiscri")
		fieldValue = 0x5;
	else if (dbgSelection == "trigger")
		fieldValue = 0x6;
	else
		fieldValue = 0x7;
	BitFieldHelper<uint32_t>::setField(regValue, 0, 3, fieldValue);

	// "read_from_0". If set, force to start the readout from column 0.
	const bool isReadFromColumn0 = agetCfg("Global")("Reg2")("isReadFromColumn0").getBoolValue();
	MSG_IF(verbose_) << "isReadFromColumn0: " << isReadFromColumn0;
	BitFieldHelper<uint32_t>::setField(regValue, 3, 1, isReadFromColumn0?0x1:0x0);

	// "test_digout". If 1, a test pattern is serialized to the output instead of the 9bit
	// address of the last read column.
	const bool isDigitalOutputTest = agetCfg("Global")("Reg2")("isDigitalOutputTest").getBoolValue();
	MSG_IF(verbose_) << "isDigitalOutputTest: " << isDigitalOutputTest;
	BitFieldHelper<uint32_t>::setField(regValue, 4, 1, isDigitalOutputTest?0x1:0x0);

	// "set_l0_when_rst". Select the analog reset level during "reset operation".
	const bool isResetLevelUndefined = agetCfg("Global")("Reg2")("isResetLevelUndefined").getBoolValue();
	MSG_IF(verbose_) << "isResetLevelUndefined";
	BitFieldHelper<uint32_t>::setField(regValue, 5, 1, isResetLevelUndefined?0x1:0x0);

	// "en_mker_rst". If 1, a "digital" marker (near gnd or vdd levels) is multiplexed to
	// the analog output during "reset operation".
	const bool isMarkerEnable = agetCfg("Global")("Reg2")("isMarkerEnable").getBoolValue();
	MSG_IF(verbose_) << "isMarkerEnable: " << isMarkerEnable;
	BitFieldHelper<uint32_t>::setField(regValue, 6, 1, isMarkerEnable?0x1:0x0);

	// "rst_lv_to1". Set the level of the digital marker (when en_mker_rst=1).
	const bool isResetLevelVdd = agetCfg("Global")("Reg2")("isResetLevelVdd").getBoolValue();
	MSG_IF(verbose_) << "isResetLevelVdd: " << isResetLevelVdd;
	BitFieldHelper<uint32_t>::setField(regValue, 7, 1, isResetLevelVdd?0x1:0x0);

	// "boost_pw". If 1, the output current of the analog block Gain-2 is increased (+20%).
	const bool isBoostPower = agetCfg("Global")("Reg2")("isBoostPower").getBoolValue();
	MSG_IF(verbose_) << "isBoostPower: " << isBoostPower;
	BitFieldHelper<uint32_t>::setField(regValue, 8, 1, isBoostPower?0x1:0x0);

	// "out_resync". If 1, the SC output data is resynchronized by a clock edge (selected by synchro_inv).
	const bool isOutputscCkSynchro = agetCfg("Global")("Reg2")("isOutputscCkSynchro").getBoolValue();
	MSG_IF(verbose_) << "isOutputscCkSynchro: " << isOutputscCkSynchro;
	BitFieldHelper<uint32_t>::setField(regValue, 9, 1, isOutputscCkSynchro?0x1:0x0);

	// "synchro_inv". select the edge for the synchronizing of the SC output data (0= rising, 1=falling).
	const bool isSynchroEdgeFalling = agetCfg("Global")("Reg2")("isSynchroEdgeFalling").getBoolValue();
	MSG_IF(verbose_) << "isSynchroEdgeFalling: " << isSynchroEdgeFalling;
	BitFieldHelper<uint32_t>::setField(regValue, 10, 1, isSynchroEdgeFalling?0x1:0x0);

	// "force_eout". If 1, inhibit the 3rd state functionality of the SC output buffer.
	const bool isSCoutputBufferInTriState = agetCfg("Global")("Reg2")("isSCoutputBufferInTriState").getBoolValue();
	MSG_IF(verbose_) << "isSCoutputBufferInTriState: " << isSCoutputBufferInTriState;
	BitFieldHelper<uint32_t>::setField(regValue, 11, 1, isSCoutputBufferInTriState?0x1:0x0);

	// "Cur_RA". These 2 bits manage the current of the SCA line buffers.
	const std::string currentBias = agetCfg("Global")("Reg2")("currentBias").getStringValue();
	MSG_IF(verbose_) << "currentBias: " << currentBias;
	if (currentBias == "current0")
		fieldValue = 0x0;
	else if (currentBias == "current1")
		fieldValue = 0x1;
	else if (currentBias == "current2")
		fieldValue = 0x2;
	else if (currentBias == "current3")
		fieldValue = 0x3;
	BitFieldHelper<uint32_t>::setField(regValue, 12, 2, fieldValue);

	// "Cur_BUF". These 2 bits manage the current of the SCA output buffers.
	const std::string currentBuffer = agetCfg("Global")("Reg2")("bufferCurrent").getStringValue();
	MSG_IF(verbose_) << "bufferCurrent: " << currentBuffer;
	if (currentBuffer == "1.503mA")
		fieldValue = 0x0;
	else if (currentBuffer == "1.914mA")
		fieldValue = 0x1;
	else if (currentBuffer == "2.7mA")
		fieldValue = 0x2;
	else if (currentBuffer == "4.870mA")
		fieldValue = 0x3;
	BitFieldHelper<uint32_t>::setField(regValue, 14, 2, fieldValue);

	// "power_down_write". If 1, put the write section in power down mode.
	const bool isPWDonWrite = agetCfg("Global")("Reg2")("isPWDonWrite").getBoolValue();
	MSG_IF(verbose_) << "isPWDonWrite: " << isPWDonWrite;
	BitFieldHelper<uint32_t>::setField(regValue, 16, 1, isPWDonWrite?0x1:0x0);

	// "power_down_read". If 1, put the read section in power down mode.
	const bool isPWDonRead = agetCfg("Global")("Reg2")("isPWDonRead").getBoolValue();
	MSG_IF(verbose_) << "isPWDonRead: " << isPWDonRead;
	BitFieldHelper<uint32_t>::setField(regValue, 17, 1, isPWDonRead?0x1:0x0);

	// "alternate_power". If 1, set alternatively the read and write sections in power down mode.
	const bool isAlternatePower = agetCfg("Global")("Reg2")("isAlternatePower").getBoolValue();
	MSG_IF(verbose_) << "isAlternatePower: " << isAlternatePower;
	BitFieldHelper<uint32_t>::setField(regValue, 18, 1, isAlternatePower?0x1:0x0);

	// "ShortReadSeq". It control the length of the "reset level" in the SCA readout phase.
	const bool is1bitShortResetReadSeq = agetCfg("Global")("Reg2")("is1bitShortResetReadSeq").getBoolValue();
	MSG_IF(verbose_) << "is1bitShortResetReadSeq: " << is1bitShortResetReadSeq;
	BitFieldHelper<uint32_t>::setField(regValue, 19, 1, is1bitShortResetReadSeq?0x1:0x0);

	// "Dis_Multiplicity_Out". If 1, disable the trigger output.
	const bool isTriggerOutputDisable = agetCfg("Global")("Reg2")("isTriggerOutputDisable").getBoolValue();
	MSG_IF(verbose_) << "isTriggerOutputDisable: " << isTriggerOutputDisable;
	BitFieldHelper<uint32_t>::setField(regValue, 20, 1, isTriggerOutputDisable?0x1:0x0);

	// "Autoreset Bank". Select the reset mode of the SCA pointer.
	const bool isRisingedgeWriteAutoresetBank = agetCfg("Global")("Reg2")("isRisingedgeWriteAutoresetBank").getBoolValue();
	MSG_IF(verbose_) << "isRisingedgeWriteAutoresetBank: " << isRisingedgeWriteAutoresetBank;
	BitFieldHelper<uint32_t>::setField(regValue, 21, 1, isRisingedgeWriteAutoresetBank?0x1:0x0);

	// "En_trigg_lvds". if 1, enable the lvds trigger signal on the output pads 85 & 86.
	const bool islvdsTriggerOutput = agetCfg("Global")("Reg2")("islvdsTriggerOutput").getBoolValue();
	MSG_IF(verbose_) << "islvdsTriggerOutput: " << islvdsTriggerOutput;
	BitFieldHelper<uint32_t>::setField(regValue, 22, 1, islvdsTriggerOutput?0x1:0x0);

	// Bit "Input_dynamic_range". See AGET Data Sheet, January 2013
	// The input dynamic range of the discriminator can be fixed to 5% (1) or 17.5% (0) of the input dynamic range of the analog channel.
	bool useInputDynamicRange5pct = false; // Default value of 0 for compatibility with AGET 2.1 where this bit does no exist.
	try
	{
		useInputDynamicRange5pct = agetCfg("Global")("Reg2")("inputDynamicRange_5pct").getBoolValue();
	}
	catch (const CCfg::CConfig::NotFound &)
	{
		// keep default value
	}
	const uint32_t bitValue = useInputDynamicRange5pct ? 0x1 : 0x0;
	MSG_IF(verbose_) << "Input_dynamic_range: " << bitValue;
	BitFieldHelper<uint32_t>::setField(regValue, 24, 1, bitValue);

	// Bits 25-26 "trigger_LVDS_level". Specific to ASTRE chip, added Nov 2016.
	// Controls multiplicity level of LVDS trigger, if bit 22 is set.
	uint32_t pow_lvds_level = (uint32_t) agetCfg("Global")("Reg2")("triggerLvdsThreshold").getIntValue(1u);
	uint32_t lvds_level = logtwo(pow_lvds_level);
	MSG_IF(verbose_) << "LVDS trigger threshold: " << lvds_level;
	BitFieldHelper<uint32_t>::setField(regValue, 25, 2, lvds_level);

	// Bits 27-28 "unity_multiplicity_level". Specific to ASTRE chip, added Nov 2016.
	// Controls unity multiplicity (e.g. per channel).
	// AGET behavior for 0; for 1, 32 channels cover all range; for 2, 16 channels; or 3, 8 channels.
	uint32_t pow_unit_mult = (uint32_t) agetCfg("Global")("Reg2")("multiplicityLevel").getIntValue(1u);
	uint32_t unit_mult = logtwo(pow_unit_mult);
	MSG_IF(verbose_) << "Multiplicity unity level: " << unit_mult;
	BitFieldHelper<uint32_t>::setField(regValue, 27, 2, unit_mult);

	return regValue;
}
//_________________________________________________________________________________________________
/**
 * Parses a configuration file and computes the desired value for AGET 16-bit register 12.
 */
uint16_t CoBoConfigParser::agetReg12Value(CCfg::CConfig & cfg, const size_t & asadId,
		const size_t & agetId)
{
	const uint16_t regValue = agetConfig(cfg, asadId, agetId)("Global")("ReadoutPointerOffset").getIntValue();
	MSG_IF(verbose_) << "Value used for register 'reg12': " << std::showbase << std::hex << regValue;
	return regValue;
}
//_________________________________________________________________________________________________
/**
 * Parses a configuration file and computes the desired value for an AGET register.
 * @param regName Name of AGET register to fill.
 * @param cfg Configuration for the board.
 * @param asadId AsAd index.
 * @param agetId AGET index.
 * @return Value of the specified register corresponding to configuration.
 */
uint64_t CoBoConfigParser::agetRegValue(CCfg::CConfig cfg, const size_t & asadId, const size_t & agetId, const std::string & regName, const uint16_t & agetVersion)
{
	if ("reg1" == regName)
	{
		return agetReg1Value(cfg, asadId, agetId, agetVersion);
	}
	else if ("reg2" == regName)
	{
		return agetReg2Value(cfg, asadId, agetId);
	}
	else if ("reg12" == regName)
	{
		return agetReg12Value(cfg, asadId, agetId);
	}

	// General case
	uint64_t regValue = 0;

	CCfg::CConfig::Iterator itChan = agetConfig(cfg, asadId, agetId).iterate("channel");

	while (itChan.hasNext())
	{
		itChan.next();
		if (CCfg::DEFAULT_OBJ_INDEX == itChan->getIndex())
			continue;
		if ((*itChan)("isActive").getBoolValue())
		{
			size_t channelId = boost::lexical_cast<size_t>(itChan->getIndex());
			uint64_t bitMask = agetChanBitmask(regName, channelId, itChan);
			regValue |= bitMask;
		}
	}

	return regValue;
}
//_________________________________________________________________________________________________
size_t CoBoConfigParser::agetChannelLsbThresholdOffset(const std::string & regName, const size_t & channelId)
{
	size_t offset = 64;
	if (regName == "reg8msb")
	{
		// Builds mask for setting channel individual threshold bits in AGET register 8.
		// When numbering channels from 1 to 64 without FPN, the 64 MS bits of register 8 deal with channels 16 to 1
		// When numbering channels from 0 to 67 with    FPN, the 64 MS bits of register 8 deal with channels [16-12]U[10-0]
		size_t physicalChannelId = 1 + agetPhysicalChannelId(channelId);
		if (physicalChannelId <= 16 and physicalChannelId >= 1)
		{
			offset = 64 - 4 * physicalChannelId;
		}
	} else if (regName == "reg8")
	{
		// Builds mask for setting channel individual threshold bits in AGET register 8.
		// When numbering channels from 1 to 64 without FPN, the 64 LS bits of register 8 deal with channels 32 to 17
		// When numbering channels from 0 to 67 with    FPN, the 64 LS bits of register 8 deal with channels [33-23]U[21-17]
		size_t physicalChannelId = 1 + agetPhysicalChannelId(channelId);
		if (physicalChannelId <= 32 and physicalChannelId >= 17)
		{
			offset = 64 - 4 * (physicalChannelId - 16);
		}
	} else if (regName == "reg9msb")
	{
		// Builds mask for setting channel individual threshold bits in AGET register 9.
		// When numbering channels from 1 to 64 without FPN, the 64 MS bits of register 9 deal with channels 49 to 64
		// When numbering channels from 0 to 67 with    FPN, the 64 MS bits of register 9 deal with channels [51-55]U[57-67]
		size_t physicalChannelId = 1 + agetPhysicalChannelId(channelId);
		if (physicalChannelId >= 49 and physicalChannelId <= 64)
		{
			offset = 4 * (physicalChannelId - 49);
		}
	} else if (regName == "reg9")
	{
		// Builds mask for setting channel individual threshold bits in AGET register 9.
		// When numbering channels from 1 to 64 without FPN, the 64 LS bits of register 9 deal with channels 33 to 48
		// When numbering channels from 0 to 67 with    FPN, the 64 LS bits of register 9 deal with channels [34-44]U[46-50]
		size_t physicalChannelId = 1 + agetPhysicalChannelId(channelId);
		if (physicalChannelId >= 33 and physicalChannelId <= 48)
		{
			offset = 4 * (physicalChannelId - 33);
		}
	}
	return offset;
}
//_________________________________________________________________________________________________
uint64_t CoBoConfigParser::agetChanBitmask(const std::string & regName, const size_t & channelId,
		CCfg::CConfig::Iterator & itChan)
{
	uint64_t bitMask = 0;

	if (regName == "reg3")
	{
		if ((*itChan)("isSelectedforTestMode").getBoolValue())
		{
			if (channelId <= 33)
			{
				int shiftIndex = 33 - channelId;
				bitMask = UINT64_C(1) << shiftIndex;
			}
		}
	} else if (regName == "reg4")
	{
		if ((*itChan)("isSelectedforTestMode").getBoolValue())
		{
			if (channelId >= 34 and channelId <= 67)
			{
				int shiftIndex = channelId - 34;
				bitMask = UINT64_C(1) << shiftIndex;
			}
		}
	} else if (regName == "reg6")
	{
		if (channelId <= 33 and channelId != 11 and channelId != 22)
		{
			int shiftIndex = 62 - 2 * channelId + 4;
			if (channelId < 11)
			{
				shiftIndex = 62 - 2 * channelId;
			} else if (channelId < 22 and channelId > 11)
			{
				shiftIndex = 62 - 2 * channelId + 2;
			}

			const std::string gain = (*itChan)("Gain").getStringValue();
			if (gain == "120fC")
			{
				bitMask = UINT64_C(0) << shiftIndex;
			} else if (gain == "240fC")
			{
				bitMask = UINT64_C(1) << shiftIndex;
			} else if (gain == "1pC")
			{
				bitMask = UINT64_C(2) << shiftIndex;
			} else if (gain == "10pC")
			{
				bitMask = UINT64_C(3) << shiftIndex;
			}
		}
	} else if (regName == "reg7")
	{
		if (channelId >= 34 and channelId <= 67 and channelId != 45 and channelId != 56)
		{
			int shiftIndex = 2 * (channelId - 34) - 4;
			if (channelId < 45)
			{
				shiftIndex = 2 * (channelId - 34);
			} else if (channelId < 56)
			{
				shiftIndex = 2 * (channelId - 34) - 2;
			}

			const std::string gain = (*itChan)("Gain").getStringValue();
			if (gain == "120fC")
			{
				bitMask = UINT64_C(0) << shiftIndex;
			} else if (gain == "240fC")
			{
				bitMask = UINT64_C(1) << shiftIndex;
			} else if (gain == "1pC")
			{
				bitMask = UINT64_C(2) << shiftIndex;
			} else if (gain == "10pC")
			{
				bitMask = UINT64_C(3) << shiftIndex;
			}
		}
	} else if (regName == "reg8" or regName == "reg8msb" or regName == "reg9" or regName == "reg9msb")
	{
		uint64_t threshold = (*itChan)("LSBThresholdValue").getIntValue();
		size_t offset = agetChannelLsbThresholdOffset(regName, channelId);
		if (offset < 64)
			bitMask = (UINT64_C(0xF) & threshold) << offset;
	} else if (regName == "reg10")
	{
		// Builds mask for setting inhibition mode of channels 0 to 33 (except for FPN channels 11 and 22) in AGET register 10.
		// See section 9.7.K of AGET Data Sheet.
		if (channelId <= 33 and channelId != 11 and channelId != 22)
		{
			// Position for channel 2 bits in register
			size_t physicalChannelIndex = channelId;
			if (channelId > 11) physicalChannelIndex--;
			if (channelId > 22) physicalChannelIndex--;
			const size_t shiftIndex = 64-2*(physicalChannelIndex + 1);

			// Value of 2 bits
			// 0x0 = no inhibition
			// 0x1 = channel contributes neither to multiplicity nor to hit register
			// 0x2 = channel does not contribute to multiplicity
			// 0x3 is equivalent to 0x1
			std::string inhibitionMode = (*itChan)("TriggerInhibition").getStringValue();
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
			BitFieldHelper<uint64_t>::setField(bitMask, shiftIndex, 2, modeValue);
		}
	} else if (regName == "reg11")
	{
		// Builds mask for setting inhibition mode of channels 34 to 67 (except for FPN channels 45 and 56) in AGET register 11.
		// See section 9.7.L of AGET Data Sheet.
		if (channelId >= 34 and channelId <= 67 and channelId != 45 and channelId != 56)
		{
			// Position for channel 2 bits in register
			size_t physicalChannelIndex = channelId;
			if (channelId > 45) physicalChannelIndex--;
			if (channelId > 56) physicalChannelIndex--;
			const size_t shiftIndex = 2*(physicalChannelIndex-34);

			// Value of 2 bits
			// 0x0 = no inhibition
			// 0x1 = channel contributes neither to multiplicity nor to hit register
			// 0x2 = channel does not contribute to multiplicity
			// 0x3 is equivalent to 0x1
			std::string inhibitionMode = (*itChan)("TriggerInhibition").getStringValue();
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
			BitFieldHelper<uint64_t>::setField(bitMask, shiftIndex, 2, modeValue);
		}
	}
	else
	{
		throw std::string("Register \'") + regName + "\' does not exist";
	}

	return bitMask;
}
//_________________________________________________________________________________________________
/**
 * Reads the board configuration and creates a mask to set the registers indicating which AGET channels to read.
 * @param cfg Configuration for the board.
 * @param asadId Index of the AsAd board.
 * @param agetId Index of the AGET.
 * @param maskName Mask name: "only_if_hit" of "always"
 * @param minChanId
 * @param maxChanId
 */
uint32_t CoBoConfigParser::agetReadoutMask(CCfg::CConfig cfg, const size_t & asadId, const size_t & agetId,
		const std::string & maskName, const size_t minChanId, const size_t & maxChanId)
{
	uint32_t regValue = 0;
	// Loop over configurations of concerned channels
	for (size_t chanId=minChanId; chanId <= maxChanId; ++chanId)
	{
		CCfg::CConfig chanCfg(agetConfig(cfg, asadId, agetId)("channel", chanId));
		if (chanCfg("isActive").getBoolValue())
		{
			const std::string chanParam = chanCfg("Reading").getStringValue();
			if (maskName == chanParam)
			{
				BitFieldHelper<uint32_t>::setBit(regValue, maxChanId - chanId, true);
			}
		}
	}
	return regValue;
}
//_________________________________________________________________________________________________
/**
 * Extracts from the configuration the value of a parameter of the circular buffer.
 */
uint32_t CoBoConfigParser::getCircularBufferParameter(CCfg::CConfig cfg, const std::string & parName)
{
	return cfg("CircularBuffer")(parName).getHexValue();
}
//_________________________________________________________________________________________________
/**
 * Extracts from the configuration the value of a parameter of the circular buffer of a specific AsAd board.
 */
uint32_t CoBoConfigParser::getAsAdCircularBufferParameter(CCfg::CConfig cfg, const size_t & asadId, const std::string & parName)
{
	return asadConfig(cfg, asadId)("CircularBuffer")(parName).getHexValue();
}
//_________________________________________________________________________________________________
/**
 * Returns 4-bit integer defining the function of the LEMO connector on the CoBo front panel.
 *
 * The lemoMux pin configures each of the LEMO I/O connectors on the front panel of the CoBo.
 * Valid options for each connector are:
 * - 0x0: trigger pulse output
 * - 0x1: busy output
 * - 0x2: pll_holdover output
 * - 0x3: pll_ld output
 * - 0x4: pll_clkin0 output
 * - 0x5: pll_clkin1 output
 * - 0x6: pll_sync output
 * - 0x8: external trigger input
 * - 0xC: pll_clkin0 input
 * - 0xD: pll_clkin1 input
 * - 0xE: pll_sync input
 * - 0xF: external timestamp reset input
 * - others: unused
 *
 * The firmware defaults to value 0x9 for all LEMOs, so they will not function unless configured by slow control.
 * It is ok for more than one LEMO to share the same function.
 * If multiple LEMOs are configured to the same input signal, the firmware will OR the signals.
 *
 * @param cfg CoBo board configuration.
 * @param lemoIndex Index of the LEMO connector (0 to 3).
 */
uint8_t CoBoConfigParser::getLemoConnectorMode(CCfg::CConfig cfg, const size_t & lemoIndex)
{
	uint8_t mode = 0x9; // default firmware value
	// List of functions of LEMO connector
	static const std::map< std::string, uint8_t > lemoModes
	 	 = boost::assign::map_list_of
	 	("trigger_pulse_output", 0x0)
		("busy_output", 0x1)
		("pll_holdover_output", 0x2)
		("pll_ld_output", 0x3)
		("pll_clkin0_output", 0x4)
		("pll_clkin1_output", 0x5)
		("pll_sync_output", 0x6)
		("external_trigger_input", 0x8)
		("none", 0x9)
		("pll_clkin0_input", 0xC)
		("pll_clkin1_input", 0xD)
		("pll_sync_input", 0xE)
		("external_timestamp_reset_input", 0xF);
	try
	{
		// Get string chosen by user
		std::string modeString = cfg("LEMO")("L", lemoIndex)("function").getStringValue();
		// Check it corresponds to a valid mode
		std::map< std::string, uint8_t >::const_iterator modeIter;
		modeIter = lemoModes.find(modeString);
		if (lemoModes.end() != modeIter)
		{
			mode = modeIter->second;
		}
		else
		{
			LOG_ERROR() << "'" << modeString << "' is not a valid LEMO connector mode!";
		}
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		// Added for M-CoBo firmware based on N. U.'s email of May 14th, 2014 (firmwareDate = 0x20140514, firmwareRelease = 7)
		LOG_WARN() << e.what();
	}

	return mode;
}
//_________________________________________________________________________________________________
/**
 * Returns, for each AsAd board, the identifier of the chip to use for serial delay calibration.
 * @return 8-bit field with 2bits per AsAd board containing the chip identifier between 0 and 3.
 */
uint8_t CoBoConfigParser::getCalibrationChips(CCfg::CConfig cfg)
{
	uint8_t userMask = 0;
	for (size_t asadId=0; asadId <4; asadId++)
	{
		try
		{
			uint8_t chipId = asadConfig(cfg, asadId)("Control")("calibrationAgetId").getIntValue();
			BitFieldHelper< uint8_t >::setField(userMask, asadId*2, 2u, chipId);
		}
		catch (const CCfg::CConfig::NotFound & e)
		{
			//LOG_WARN() << e.what();
		}
	}
	return userMask;
}
//_________________________________________________________________________________________________
/**
 * Returns default value of PLL register, modified with parameters parsed from the configuration file.
 * The default values are taken from PLL configuration in email from F. A.-N. on May 16th, 2014
 * modified according to email from F. A.-N. on March 17th, 2015.
 * @param cfg CoBo board configuration.
 * @param regIndex Index of the PLL register.
 */
uint32_t CoBoConfigParser::getPllRegValue(CCfg::CConfig cfg, const size_t & regIndex)
{
	uint32_t regValue = 0;
	try
	{
		switch (regIndex)
		{
		case 10u:
			// R11: VCO_DIV [10:8]
			regValue = 0x9140410A;
			utl::BitFieldHelper< uint32_t >::setField(regValue, 8u, 3u, cfg("Module")("PLL")("VCO_DIV").getIntValue());
			break;
		case 11u:
			// R11: MODE [31:27], EN_SYNC [26], SYNC_POL_INV [16], SYNC_EN_AUTO [15], SYNC_TYPE [14:12]
			regValue = 0x0400000B;
			utl::BitFieldHelper< uint32_t >::setField(regValue, 27u, 5u, cfg("Module")("PLL")("MODE").getIntValue());
			utl::BitFieldHelper< uint32_t >::setBit(regValue, 26u, cfg("Module")("PLL")("EN_SYNC").getBoolValue());
			utl::BitFieldHelper< uint32_t >::setBit(regValue, 16u, cfg("Module")("PLL")("SYNC_POL_INV").getBoolValue());
			utl::BitFieldHelper< uint32_t >::setBit(regValue, 15u, cfg("Module")("PLL")("SYNC_EN_AUTO").getBoolValue());
			utl::BitFieldHelper< uint32_t >::setField(regValue, 12u, 3u, cfg("Module")("PLL")("SYNC_TYPE").getIntValue());
			break;
		case 12u:
			// R12: LD_MUX [31:27], LD_TYPE [26:24], SYNC_PLL2_DLD [23], SYNC_PLL1_DLD [22]
			regValue = 0x0B0C006C;
			utl::BitFieldHelper< uint32_t >::setField(regValue, 27u, 5u, cfg("Module")("PLL")("LD_MUX").getIntValue());
			utl::BitFieldHelper< uint32_t >::setField(regValue, 24u, 3u, cfg("Module")("PLL")("LD_TYPE").getIntValue());
			try
			{
				// Added after email from F. A.-N. on March 16th, 2015
				utl::BitFieldHelper< uint32_t >::setBit(regValue, 23u, cfg("Module")("PLL")("SYNC_PLL2_DLD").getBoolValue());
				utl::BitFieldHelper< uint32_t >::setBit(regValue, 22u, cfg("Module")("PLL")("SYNC_PLL1_DLD").getBoolValue());
			}
			catch (const CCfg::CConfig::NotFound & e)
			{
				LOG_WARN() << e.what();
			}
			break;
		case 13u:
			// R13: HOLDOVER_MUX [31:27], HOLDOVER_TYPE [26:24], Status_CLKin1_MUX [22:20], Status_CLKin0_TYPE [18:16], Status_CLKin0_MUX [14:12], CLKin_Select_MODE [11:9], EN_CLKin1 [6], EN_CLKin0 [5]
			regValue = 0x1B12902D;
			utl::BitFieldHelper< uint32_t >::setField(regValue, 27u, 5u, cfg("Module")("PLL")("HOLDOVER_MUX").getIntValue());
			utl::BitFieldHelper< uint32_t >::setField(regValue, 24u, 3u, cfg("Module")("PLL")("HOLDOVER_TYPE").getIntValue());
			utl::BitFieldHelper< uint32_t >::setField(regValue, 20u, 3u, cfg("Module")("PLL")("Status_CLKin1_MUX").getIntValue());
			utl::BitFieldHelper< uint32_t >::setField(regValue, 16u, 3u, cfg("Module")("PLL")("Status_CLKin0_TYPE").getIntValue());
			utl::BitFieldHelper< uint32_t >::setField(regValue, 12u, 3u, cfg("Module")("PLL")("Status_CLKin0_MUX").getIntValue());
			utl::BitFieldHelper< uint32_t >::setField(regValue,  9u, 3u, cfg("Module")("PLL")("CLKin_SELECT_MODE").getIntValue());
			try
			{
				// Added after email from F. A.-N. on March 16th, 2015
				utl::BitFieldHelper< uint32_t >::setBit(regValue, 6u, cfg("Module")("PLL")("EN_CLKin1").getBoolValue());
				utl::BitFieldHelper< uint32_t >::setBit(regValue, 5u, cfg("Module")("PLL")("EN_CLKin0").getBoolValue());
			}
			catch (const CCfg::CConfig::NotFound & e)
			{
				LOG_WARN() << e.what();
			}
			break;
		case 14u:
			// R14: Status_CLKin1_TYPE [26:24]
			regValue = 0x1200000E;
			utl::BitFieldHelper< uint32_t >::setField(regValue, 24u, 3u, cfg("Module")("PLL")("Status_CLKin1_TYPE").getIntValue());
			break;
		default:
			break;
		}
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	return regValue;
}
//_________________________________________________________________________________________________
/**
 * Returns (initial) value of the difference between output voltage and reference voltage of the AsAd generator.
 *
 * @param cfg CoBo board configuration.
 * @param asadId AsAd board identifier.
 * @return Value of the difference between output voltage and reference voltage of the AsAd generator [mV].
 */
int16_t CoBoConfigParser::asadPulserVoltage_mV(CCfg::CConfig cfg, const size_t & asadId)
{
	int16_t amplitude = 0;
	try
	{
		amplitude = (int16_t) cfg("Generator")("amplitudeStart").getValueAsReal();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		LOG_WARN() << "Trying former configuration file format with AsAd no. 0";
		amplitude = asadConfig(cfg, asadId)("Generator")("Amplitude", 0).getValueAsReal();
	}
	return amplitude;
}
//________________________________________________________________________________________________
/**
 * Extracts value for AsAd generator capacitor bit RNG (aka Rg) of register TCM0 of AsAd TCM device.
 *
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
 * @param cfg CoBo board configuration
 * @param asadId AsAd board identifier
 */
bool CoBoConfigParser::asadPulser11pFCapacitorSelected(CCfg::CConfig & cfg, const size_t & asadId)
{
	// Get user parameter value
	std::string capacitanceStr = "1pF";
	try
	{
		capacitanceStr = asadConfig(cfg, asadId)("Generator")("injectValue").getStringValue();

	}
	// For a time, parameter was set at the CoBo level, identical for all four AsAd boards
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		capacitanceStr = cfg("Generator")("injectValue").getStringValue();
	}

	// Interpret string
	const bool use11pF = ("10pF" == capacitanceStr or "11pF" == capacitanceStr);
	return use11pF;
}
//________________________________________________________________________________________________
/**
 * Checks if TCM device should be configured for step voltage requests from AsA external signals or from loading the TCM2 register.
 * @see GET-QA-000-0005
 * @see GET-QA-000-0003
 *
 * @return Returns true if TRG should be set (external trigger mode), false otherwise (slow control mode).
 * @param cfg CoBo board configuration
 * @param asadId AsAd board identifier
 */
bool CoBoConfigParser::asadPulserTriggerMode(CCfg::CConfig & cfg, const size_t & /* asadId */)
{
	// Get user parameter value
	bool externalMode = false;
	try
	{
		externalMode = cfg("Generator")("externalTriggerMode").getBoolValue();
	}
	// Parameter added on March 28th, 2014
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	return externalMode;
}
//________________________________________________________________________________________________
/**
 * Returns TCM device output step triggering delay (5 bits).
 * @see GET-QA-000-0005
 *
 * @return Delay expressed in number of WCK periods (0 to 31).
 * @param cfg CoBo board configuration
 * @param asadId AsAd board identifier
 */
uint8_t CoBoConfigParser::asadPulserTriggerDelay(CCfg::CConfig & cfg, const size_t & asadId)
{
	// Get user parameter value
	uint8_t delay_wckPeriods = 0;
	try
	{
		delay_wckPeriods = asadConfig(cfg, asadId)("Generator")("triggerDelay").getIntValue();
	}
	// Parameter added on March 28th, 2014
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	return delay_wckPeriods;
}
//________________________________________________________________________________________________
/**
 * Checks mode to select for the AsAd generator in the TCM0 register.
 *  - Mode "test": a step voltage goes directly to the AGET test input.
 *  - Mode "calibration": a current pulse goes through a capacitor to the AGET calibration input.
 *
 * @see Section 1.4 of GET-AS-002-0005: "AsAd Test and Calibration Management".
 *
 * @return Returns true for calibration mode, false for test mode.
 * @param cfg CoBo board configuration
 * @param asadId AsAd board identifier
 */
bool CoBoConfigParser::asadPulserCalibMode(CCfg::CConfig & cfg, const size_t & asadId)
{
	bool isCalibMode = false;
	try
	{
		isCalibMode = asadConfig(cfg, asadId)("Generator")("isCalibMode").getBoolValue();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		isCalibMode = cfg("Generator")("isCalibMode").getBoolValue();
	}
	return isCalibMode;
}
//________________________________________________________________________________________________
/**
 * Returns the difference between default output voltage 1.102 mV.
 *
 * @param cfg CoBo board configuration.
 * @return Value of the difference between default output voltage and reference voltage of the AsAd generator [mV].
 */
float CoBoConfigParser::asadPulserDefaultVoltage_mV(CCfg::CConfig & cfg)
{
	float defaultVoltage = 0;
	try
	{
		defaultVoltage = cfg("Generator")("defaultVoltage").getValueAsReal();
	}
	catch (CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
	}
	return defaultVoltage;
}
//________________________________________________________________________________________________
/**
 * Parse settings for the AsAd inspection manager.
 */
void CoBoConfigParser::getAsadInspectionLinesParameters(CCfg::CConfig & cfg, const size_t & asadId,
		std::string & deviceName, std::string & isp1Name, std::string & isp2Name)
{
	asadConfig(cfg, asadId)("Generator")("isCalibMode").getBoolValue();
	try
	{
		deviceName = asadConfig(cfg, asadId)("InspectionLines")("deviceToInspect").getStringValue();
		isp1Name = asadConfig(cfg, asadId)("InspectionLines")(deviceName)("ISP1").getStringValue();
		isp2Name = asadConfig(cfg, asadId)("InspectionLines")(deviceName)("ISP2").getStringValue();
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		LOG_WARN() << e.what();
		deviceName = "ADC";
		isp1Name = "SPI_CS";
		isp2Name = "CAS";
	}
}
//_________________________________________________________________________________________________
/**
 * Returns identifier of single channel item used for AsAd monitoring in single channel mode.
 * 000: VDD1; 001: Internal temperature; 010: External temperature; 100: VDD; 101: IDD
 */
uint8_t CoBoConfigParser::getAsAdMonitoringSingleChannel(CCfg::CConfig & cfg, const size_t & asadIdx)
{
	const std::string itemString = asadConfig(cfg, asadIdx)("Monitoring")("singleChannelItem").getStringValue();
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
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
