/**
 * @file CoBoConfigParser.h
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

#ifndef get_rc_CoBoConfigParser_h_INCLUDED
#define get_rc_CoBoConfigParser_h_INCLUDED

#include "CCfg/CConfig.h"
#include <boost/cstdint.hpp>
#include <string>

namespace get
{
namespace rc
{

class CoBoConfigParser
{
public:
	CoBoConfigParser();
	static bool isAsadActive(const CCfg::CConfig& cfg, size_t asadId);
	static bool isSCRAutoDelayEnabled(CCfg::CConfig cfg, size_t asadId);
	static CCfg::CConfig asadConfig(CCfg::CConfig cfg, const size_t & asadId);
	static CCfg::CConfig agetConfig(CCfg::CConfig cfg, const size_t & asadId, const size_t & agetId);
	static size_t agetPhysicalChannelId(const size_t channelId);
	static uint16_t peakingTimeIndex(size_t peakingTime, uint16_t agetVersion);
	static uint32_t agetReg1Value(CCfg::CConfig & cfg, const size_t & asadId, const size_t & agetId, const uint16_t & agetVersion);
	static uint32_t agetReg2Value(CCfg::CConfig & cfg, const size_t & asadId, const size_t & agetId);
	static uint16_t agetReg12Value(CCfg::CConfig & cfg, const size_t & asadId, const size_t & agetId);
	static uint64_t agetRegValue(CCfg::CConfig cfg, const size_t & asadId, const size_t & agetId, const std::string & regName, const uint16_t & agetVersion);
	static size_t agetChannelLsbThresholdOffset(const std::string & regName, const size_t & channelId);
	static uint64_t agetChanBitmask(const std::string & regName, const size_t & channelId, CCfg::CConfig::Iterator & itChan);
	static uint32_t agetReadoutMask(CCfg::CConfig cfg, const size_t & asadId, const size_t & agetId, const std::string & maskName, const size_t minChanId, const size_t & maxChanId);
	uint32_t getAsAdCircularBufferParameter(CCfg::CConfig cfg, const size_t & asadId, const std::string & parName);
	static uint32_t getCircularBufferParameter(CCfg::CConfig cfg, const std::string & parName);
	static uint8_t getLemoConnectorMode(CCfg::CConfig cfg, const size_t & lemoIndex);
	static uint8_t getCalibrationChips(CCfg::CConfig cfg);
	static uint32_t getPllRegValue(CCfg::CConfig cfg, const size_t & regIndex);
	/// @name AsAd TCM device
	/// @{
	static float asadPulserDefaultVoltage_mV(CCfg::CConfig & cfg);
	static int16_t asadPulserVoltage_mV(CCfg::CConfig cfg, const size_t & asadId=0);
	static bool asadPulser11pFCapacitorSelected(CCfg::CConfig & cfg, const size_t & asadId);
	static bool asadPulserCalibMode(CCfg::CConfig & cfg, const size_t & asadId);
	static bool asadPulserTriggerMode(CCfg::CConfig & cfg, const size_t & asadId);
	static uint8_t asadPulserTriggerDelay(CCfg::CConfig & cfg, const size_t & asadId);
	/// @}
	static void getAsadInspectionLinesParameters(CCfg::CConfig & cfg, const size_t & asadId, std::string & deviceName, std::string & isp1Name, std::string & isp2Name);
	static uint8_t getAsAdMonitoringSingleChannel(CCfg::CConfig & cfg, const size_t & asadIdx);
	static void setVerbose(bool enable);
private:
	static bool verbose_;
};

} /* namespace rc */
} /* namespace get */
#endif /* get_rc_CoBoConfigParser_h_INCLUDED */
