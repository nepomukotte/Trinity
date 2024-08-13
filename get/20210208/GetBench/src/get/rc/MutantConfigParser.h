/**
 * @file MutantConfigParser.h
 * @date Feb 26, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: MutantConfigParser.h 1794 2021-02-03 11:41:04Z psizun $
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

#ifndef get_rc_MutantConfigParser_h_INCLUDED
#define get_rc_MutantConfigParser_h_INCLUDED

#include <CCfg/CConfig.h>
#include <map>
#include <boost/cstdint.hpp>

namespace get
{
namespace mutant
{
	typedef enum {CsDisabled, CENTRUM, Scalers} CentrumScalerMode;
}
namespace rc
{
//_________________________________________________________________________________________________
/**
 * A class to parse a Mutant module configuration file and store its parameters.
 *
 * @see Mutant Slow Control Registers (v1.4)
 */
class MutantConfigParser
{
public:
	MutantConfigParser();
	void parseConfig(CCfg::CConfig cfg);
	bool isSlave; 			///< Whether module is a master or a slave
	bool isSlavePortA; 		///< Whether slave is connected to port A or B
	size_t slaveId; 		///< Slave identifier
	bool vhdciPortAEnabled; ///< Whether VHDCI port A is enabled
	bool vhdciPortBEnabled; ///< Whether VHDCI port B is enabled

	uint64_t timestamp_10ns; ///< 48-bit value to preload in the timestamp counter, in units of 10 ns
	uint32_t eventNumber;    ///< 32-bit value to preload in the event counter
	bool preLoadTimestamp; ///< Whether to pre-load timestamp
	bool preLoadEventNumber; ///< Whether to pre-load a value into event counter

	bool l0l1; ///< Value to give to the L0/L1 bit

	bool level0Enabled; 	///< Whether to enable trigger level 0
	uint16_t level0GateDelay_10ns; ///< Delay value of the Level 0 gate, in units of 10 ns
	uint16_t level0GateDuration_10ns; ///< Duration value of the Level 0 gate, in units of 10 ns

	bool level1Enabled; 	///< Whether to enable trigger level 1
	uint8_t level1Selection; ///< Value of L1SEL bits
	uint16_t level1aLowThreshold; ///< Level 1 low threshold multiplicity trigger value A
	uint16_t level1aHighThreshold; ///< Level 1 high threshold multiplicity trigger value A
	uint16_t level1bLowThreshold; ///< Level 1 low threshold multiplicity trigger value B (2p decay mode)
	uint16_t level1bHighThreshold; ///< Level 1 high threshold multiplicity trigger value B (2p decay mode)
	uint16_t level1aLowNumBuckets; ///< Level 1A low threshold number of buckets
	uint16_t level1aHighNumBuckets; ///< Level 1A high threshold number of buckets
	uint16_t level1bLowNumBuckets; ///< Level 1B low threshold number of buckets (2p decay mode)
	uint16_t level1bHighNumBuckets; ///< Level 1B high threshold number of buckets (2p decay mode)
	uint16_t level1bTimeout_2560ns; ///< Timeout value in case of no L1B trigger, in units of 10 ns
	uint16_t level1PostTrigDelay_10ns; ///< Time added to a valid L1 trigger in order to delay the STOP sampling, in units of 10 ns

	bool level2Enabled; 	///< Whether to enable trigger level 2
	bool hitPatternModified; ///< Tells the trigger/timestamper state machine that a mask pattern will be provided to CoBos in case of L2 OK

	bool tpdModeEnabled; 	///< Whether the two proton decay mode is enabled
	bool halfEventReadoutEnabled; ///< Whether to readout half-events after timeout in 2p mode
	mutant::CentrumScalerMode csMode; ///< Whether to disable the C/S interface, or enable scaler channels or enable coupling with CENTRUM.
	bool testModeEnabled; 	///< Whether to enable the test mode.
	bool readoutModeEnabled; 	///< Whether to enable the test mode.
	bool clearScalersEnabled; ///< Whether to clear scaler counters at run start
	bool clearCountersEnabled; ///< Whether to clear event counters at run start
	std::map< size_t, bool > enableScalers; ///< Whether to enable each of the scalers (1 to 5)

	uint8_t vcxoControl;  ///< 8-bit value to assign to the VCXO frequency control; Default should be 0x80.

	uint8_t sourceId;		///< Byte included in data frames in order to identify the detection system of which Mutant is a part of.
	uint8_t frameType; 		///< Type of frames

	// AMC Tag Only
	uint16_t amcTagOnly;  ///< 12-bit value indicating which AMC slots should only receive sync info

	// User register
	uint16_t userReg;  ///< 16-bit user defined value

	// Global Master Clock
	bool externalGmcClockEnabled;
	uint8_t externalGmcClockSelection;

	// Inspection
	uint8_t frontPanelInsp1;
	uint8_t frontPanelInsp2;
	uint8_t rearPanelInsp1;
	uint8_t rearPanelInsp2;
};
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
#endif /* get_rc_MutantConfigParser_h_INCLUDED */
