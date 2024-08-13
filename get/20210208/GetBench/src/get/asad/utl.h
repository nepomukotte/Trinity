/**
 * @file utl.h
 * @date Dec 20, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: utl.h 1724 2017-11-14 08:45:16Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the GetBench software project.
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

#ifndef get_asad_utl_h_INCLUDED
#define get_asad_utl_h_INCLUDED

#include <CCfg/CConfig.h>
#include <string>
#include <boost/cstdint.hpp>

namespace get
{
namespace asad
{
namespace utl
{
//_________________________________________________________________________________________________
struct Lcm1Config
{
	uint64_t findiv; ///< 7 bits
	uint64_t fbdiv; ///< 7 bits
	uint64_t oadiv; ///< 5 bits
	uint64_t oamux; ///< 3 bits
	uint64_t fbsel; ///< 2 bits
	uint64_t fbdly; ///< 5 bits
	uint64_t xdlysel; ///< 1 bit
	uint64_t dlygla; ///< 5 bits
	uint64_t vcosel; ///< 2 bits
};
//_________________________________________________________________________________________________
struct Lcm2Config
{
	uint64_t findiv; ///< 7 bits
	uint64_t fbdiv; ///< 7 bits
	uint64_t oadiv; ///< 5 bits
	uint64_t obdiv; ///< 5 bits
	uint64_t oamux; ///< 2 bits
	uint64_t obmux; ///< 2 bits
	uint64_t fbsel; ///< 2 bits
	uint64_t fbdly; ///< 5 bits
	uint64_t xdlysel; ///< 1 bit
	uint64_t dlygla; ///< 5 bits
	uint64_t dlyglb; ///< 5 bits
};
//_________________________________________________________________________________________________
uint64_t buildIspmRegisterValue(const std::string & moduleName, const std::string & signalNameISP1, const std::string & signalNameISP2);

uint64_t buildLcm1RegisterValue(const float & freqCKW_MHz, const float & clockDivider = 1.0);
void parseAsAdLcm1DebugConfig(CCfg::CConfig cfg, asad::utl::Lcm1Config & lcm1);
void patchLcm1RegisterValue(uint64_t & lcm1Value, const Lcm1Config & config);

uint64_t buildLcm2RegisterValue(const float & freqCKR_MHz, const int32_t & adcClockPhaseShift, const int32_t & adcClockDelay);
void parseAsAdLcm2DebugConfig(CCfg::CConfig cfg, asad::utl::Lcm2Config & lcm2);
void patchLcm2RegisterValue(uint64_t & lcm2Value, const Lcm2Config & config);

std::string stringifyMonitoringAlertStatus(const uint8_t);
//_________________________________________________________________________________________________
} // namespace utl
} // namespace asad
} // namespace get

#endif /* get_asad_utl_h_INCLUDED */
