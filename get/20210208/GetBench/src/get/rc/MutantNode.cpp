/**
 * @file MutantNode.cpp
 * @date Sep 18, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: MutantNode.cpp 1799 2021-02-08 15:42:53Z psizun $
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

#include "mdaq/rc/HardwareNode.h"
#include "mdaq/hw/RemoteRegisterProcessor.h"
#include "utl/BitFieldHelper.hpp"
#include "utl/Logging.h"
#include ICE_HEADER(Communicator.h)
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <rc/MutantNode.h>
using mdaq::utl::CmdException;
using boost::lexical_cast;
using namespace get::mutant;
using namespace mdaq::hw;

#define CATCH_MUTANT_EXCEPTIONS()						\
	catch (const CmdException & e)	             		\
	{													\
		std::string msg = std::string("Error configuring MuTanT module. ") + e.reason; \
		throw CmdException(msg);						\
	}													\
	catch (const std::exception & e)					\
	{													\
		std::string msg = std::string("Error configuring MuTanT module. ") + e.what(); \
		throw CmdException(msg);						\
	}

namespace get
{
namespace rc
{

const std::string MutantNode::NODE_TYPE("Mutant");

//_________________________________________________________________________________________________
/**
 * Constructs a manager for a MuTanT module.
 */
MutantNode::MutantNode(const std::string& moduleId_, mdaq::rc::HardwareNode* node, get::rc::BackEnd* backend)
	: Node(NodeId(NODE_TYPE, moduleId_), backend, node), mutantModuleCount_(0), fMaskGmc_(0), fMaskWsca_(0)
{
	if (node)
	{
		NODE_DEBUG() << "Creating manager for MuTanT module \'" << id().toString() << '\'';
	}
}

//_________________________________________________________________________________________________
bool MutantNode::isSlaveModule() const
{
	return config.isSlave;
}
//_________________________________________________________________________________________________
// BR 13 may 2015
bool MutantNode::isReadoutEnabled() const
{
	return config.readoutModeEnabled;
}
//_________________________________________________________________________________________________
// BR 13 may 2015
bool MutantNode::isTriggerEnabled() const
{
	if (config.level0Enabled || config.level1Enabled || config.level2Enabled)
		return true;
	else
		return false;
}
//_________________________________________________________________________________________________
void MutantNode::setModuleCount(size_t n)
{
	mutantModuleCount_ = n;
}
//_________________________________________________________________________________________________
DataSenderId MutantNode::dataSenderId() const
{
	return DataSenderId(id().toString());
}
//_________________________________________________________________________________________________
DataSenderSet MutantNode::dataSenders() const
{
	DataSenderSet senders;
	senders.insert(dataSenderId());
	return senders;
}
//_________________________________________________________________________________________________
/**
 * Prepares MuTanT module.
 *  - Checks power supply of boards A and B.
 *  - Checks module identifier.
 *  - Ensures module is not running.
 */
void MutantNode::impl_prepare(CCfg::CConfig nodeCfg)
{
	NODE_INFO() << "Preparing MuTanT module '" << id().toString() << "'...";
	try
	{
		config.parseConfig(nodeCfg);
		// Check power supply bits of boards A and B
		checkPowerSupply();
		// Check module identifier
		checkIdentifier();
		// Prior to writing to any registers, the module must have been stopped. Even though it should be the case, we make sure of it.
		//BR.
		ensureModuleNotRunning();
		// Before alignment, master and slaves should already be partially configured:
		// Tune VCXO frequency so that CoBo PLL can lock on it
		tuneVcxoFrequency(config.vcxoControl);
		// Select clock source
		selectClockSource(config.isSlave, config.isSlavePortA);
		// Set register GLOBAL_CONF
		setGlobalConfigurationRegister();
		configureTdcStatAndCardBFlags(nodeCfg);
	}
	CATCH_MUTANT_EXCEPTIONS()
}

void MutantNode::impl_unprepare()
{
}

//_________________________________________________________________________________________________
/**
 * Executes the part of the configuration of the MuTanT module for the 'configure' transition.
 *   - Configures registers of MuTanT A board.
 *   - Configures registers of MuTanT B board.
 */
void MutantNode::impl_configure(CCfg::CConfig nodeCfg)
{
	NODE_INFO() << "Configuring MuTanT module '" << id().toString() << "'...";
	config.parseConfig(nodeCfg);
	configureRegistersMutantA(nodeCfg);
	configureRegistersMutantB(nodeCfg); // BR.
	registerBeastBoards(nodeCfg);
}
//_________________________________________________________________________________________________
/**
 * Checks power supply of boards A and B.
 *
 * Checks power supply bits are set in the registers of boards A and B.
 * @see MuTanT Slow Control Registers (v1.4)
 */
void MutantNode::checkPowerSupply()
{
	NODE_DEBUG() << "Checking power supply of MuTanT A board";
	// Configure GPIO device for inputs/outputs (see email from B. R. on May 20th, 2014
	device("GPIO_SW_LEDS")->writeRegister("DIR1", 0xFCFC);
	device("GPIO_SW_LEDS")->writeRegister("DIR2", 0xFFFFFFFF);
	// Read MuTanT A register BOARD_ID_PWR_STAT
	const uint32_t BOARD_ID_PWR_STAT = device("GPIO_SW_LEDS")->readRegister("BOARD_ID_PWR_STAT");
	// Check power supply bits (24 to 30) are set
	const uint32_t powerSupplyBitsA = utl::BitFieldHelper< uint32_t >::getField(BOARD_ID_PWR_STAT, 24, 7);
	if (0x7F != powerSupplyBitsA)
	{
		std::ostringstream msg;
		msg << "Error checking MuTanT A power supply bits: register BOARD_ID_PWR_STAT = " << std::hex << std::showbase << BOARD_ID_PWR_STAT;
		throw CmdException(msg.str());
	}

	/* BR. 15/09/14 */
	NODE_DEBUG() << "Checking power supply of MuTanT B board";
	// Read MuTanT B register PWR_START_B
	const uint16_t PWR_STAT_B = device("CARD_B_EMC")->readRegister("PWR_STAT_B");
	// Check 8 power supply bits (0 to 7) are set
	const uint16_t powerSupplyBitsB = utl::BitFieldHelper< uint16_t >::getField(PWR_STAT_B, 0, 8);
	// BR. 2014/09/17 prototype returns 0xdf 
	if (0xFF != powerSupplyBitsB &&
	    0xDF != powerSupplyBitsB )
	{
		std::ostringstream msg;
		msg << "Error checking MuTanT B power supply bits: register PWR_STAT_B = " << std::hex << std::showbase << PWR_STAT_B;
		throw CmdException(msg.str());
	}/**/
}
//_________________________________________________________________________________________________
/**
 * Checks identifier MuTanT module.
 *
 * Checks board A register containing MuTanT 24-bit module identifier,
 * given by three sets of switches assembled on the bottom of the MuTanT A board.
 * @see MuTanT Slow Control Registers (v1.4)
 */
void MutantNode::checkIdentifier()
{
	// Read first 24 bits of MuTanT A register BOARD_ID_PWR_STAT
	const uint32_t BOARD_ID = device("GPIO_SW_LEDS")->readField("BOARD_ID_PWR_STAT", "mutant_id");
	NODE_DEBUG() << "Checking identifier of MuTanT module: " << std::hex << std::showbase << BOARD_ID << std::dec << std::noshowbase;
}
//_________________________________________________________________________________________________
/**
 * Ensures module is not running by un-setting the RUN bit of the TRIGGER_CTRL register of MuTanT B device CARB_B_EMC.
 */
void MutantNode::ensureModuleNotRunning()
{
	NODE_DEBUG() << "Ensuring module is not running";
	device("CARD_B_EMC")->writeField("TRIGGER_CTRL", "run", 0);
}
//_________________________________________________________________________________________________
/**
 * Returns proxy to MuTanT alignment server.
 */
mutant::AlignmentServerPrx MutantNode::alignmentServer()
{
	std::string proxyOptions;
#if ICE_INT_VERSION >= 30500
	proxyOptions = "-e 1.0";
#endif
	const std::string proxyString = hwNode().endpoint.iceProxyString("MutantAlignmentServer", "default", proxyOptions);
	Ice::ObjectPrx base = hwNode().hwServerCtrlPrx->ice_getCommunicator()->stringToProxy(proxyString);
	return mutant::AlignmentServerPrx::checkedCast(base);
}
//_________________________________________________________________________________________________
/**
 * Configures registers of MuTanT A board.
 *
 * Registers CLOCK_TREE_CONFIG1 and CLOCK_TREE_CONFIG2 are set during the alignment procedure.
 * Here, we set 2 bits:
 *  - CLOCK_TREE_CONFIG2<0> (EXT_CLK100_SEL)
 *  - CLOCK_TREE_CONFIG2<2> = CLK_PATH_SEL<2> = CLK1_INT_EXT_SEL
 *
 * @see MuTanT Slow Control Registers (v1.5)
 * @see Email from B. Raine on March 18th, 2014
 * @see Tel from G. Wittwer on March 14th, 2016: EXT_CLK100_SEL value should be tied to CLK1_INT_EXT_SEL value.
 * @see Ticket #114 and emails from TadaAki Isobe and Frederic Saillant on Dec 17, 2015
 */
void MutantNode::configureRegistersMutantA(CCfg::CConfig /* cfg */)
{
	NODE_DEBUG() << "Configuring MuTanT A registers";

	// configureTdcStatAndCardBFlags(cfg);  // now done before multi-crate alignment

	uint16_t clk_path_sel = device("GPIO_CLK_TREE")->readField("CLOCK_TREE_CONFIG2", "clk_path_sel");

	// Register CLOCK_TREE_CONFIG2 of device GPIO_CLK_TREE
	// <0>: EXT_CLK100_SEL
	// = 0: CLOCK 100 MHz from the optical transceiver input is selected
	// = 1: CLOCK 100 MHz from VHDCI or FIREWIRE connectors is selected
	utl::BitFieldHelper< uint16_t >::setBit(clk_path_sel, 0u, config.externalGmcClockEnabled);

	// <1>: CLK1 or PULSE_SEL
	// = 0: CLK1 (GMC) is distributed to the enabled AMCs (CoBo)
	// = 1: START Pulse is sent to the DSTART input of the TDC and to the enabled AMC (clock alignment)

	// <2>: CLK1_INT_EXT_SEL
	// = 0: CLK1 (GMC) is delivered by the internal VCXO
	// = 1: CLK1 (GMC) is delivered by the front panel sources (VHDCI, FIREWIRE, OPTICAL)
	utl::BitFieldHelper< uint16_t >::setBit(clk_path_sel, 2u, config.externalGmcClockEnabled);

	NODE_DEBUG() << "Setting field 'clk_path_sel' of register 'CLOCK_TREE_CONFIG2' to " << std::hex << std::showbase << clk_path_sel << std::dec;
	device("GPIO_CLK_TREE")->writeField("CLOCK_TREE_CONFIG2", "clk_path_sel", clk_path_sel);

}
//_________________________________________________________________________________________________
/**
 * Sets CLOCK_TREE_CONFIG2<3>
 * <3>: CLK2 (WSCA) or PULSE_SEL
 *      = 0: START Pulse is sent to the DSTART input of the TDC and to the enabled AMC (WSCA alignment) or SYNC pulse
 *      = 1: CLK2 is distributed to the enabled AMCs (CoBo)
 */
void MutantNode::setWscaClk2(bool enabled)
{
	uint16_t clk_path_sel = device("GPIO_CLK_TREE")->readField("CLOCK_TREE_CONFIG2", "clk_path_sel");
	utl::BitFieldHelper< uint16_t >::setBit(clk_path_sel, 3u, enabled);
	device("GPIO_CLK_TREE")->writeField("CLOCK_TREE_CONFIG2", "clk_path_sel", clk_path_sel);
}
//_________________________________________________________________________________________________
/**
 * Set CLOCK_TREE_CONFIG2<9> to start or end SYNC_PULSE procedure.
 *
 * From G. Wittwer (Dec 8, 2016):
 * Once the alignment procedure is finished, with bit CLOCK_TREE_CONFIG2<3> set back to 1, the CoBo's must not be started,
 * but bit CLOCK_TREE_CONFIG2<9> known as START_PULSE is also used in this case as SYNC_PULSE (see bottom of page 10 in Slow Control Registers Manual).
 * The PLL's of the CoBo's that receive GMC and provide different output frequencies must be aligned in phase for each frequency
 * (e.g. the 25 MHz of CoBo 'n' with the 25 MHz of CoBo 'n+1', etc...).
 * This SYNC_PULSE bit must then be set back to 0, before starting the CoBo's before starting the triggers.
 *
 * This procedure is necessary after alignment and before starting the system. On the CoBo side, during this procedure,
 * register 'mutantConfig' must be set to 1 but 'InitDone' left to 0 so that CoBo can distinguish this brief passage
 * to 1 of WSCA from a start of the acquisition.
 *
 */
void MutantNode::setSyncPulse(uint8_t sync_pulse)
{
	device("GPIO_CLK_TREE")->writeField("CLOCK_TREE_CONFIG2", "start_pulse", sync_pulse);
}
//_________________________________________________________________________________________________
/**
 * Configures register TDC_STAT_CARD_B_FLGS of device GPIO_SW_LEDS of MuTanT board A.
 *
 * - <1..0>: LED0, LED1: used for first tests of MUTANT_A board
 * - <7..4>: FLGB
 * 	  - FLGB<0>: MUTANT_B Readout FIFO Flags, = 1 -> FIFO is empty
 * 	  - FLGB<1>: MUTANT_B Readout FIFO Flags, = 1 -> FIFO is full
 * 	  - FLGB<2>: MUTANT_B Readout FIFO Flags, = 1 -> FIFO is almost full
 * 	  - FLGB<3>: not used
 * - <9..8>: GMC_SEL<1..0> = External CLK1 (GMC) source selection
 *    - 00: FireWire 1 input selected (J10)
 *    - 01: FireWire 2 input selected (J9)
 *    - 10: VHDCI source from board B
 *    - 11: Optical source (transceiver)
 * - <13..10>: TDC xxx: TDC status bits (see ACAM – TDC GPX manual)
 *    - TDC_IRQ: = 1 when measure completed
 *    - TDC_ERR: = 1 when an error occurred in the measuring process
 *    - TDC_EF1: Buffer1 Empty Flag attached to TDC STOP1, = 1 when true
 *    - TDC_EF2: Buffer2 Empty Flag attached to TDC STOP2, = 1 when true
 *
 * @see MUTanT registers mapping v1.5
 *
 */
void MutantNode::configureTdcStatAndCardBFlags(CCfg::CConfig /* cfg */)
{
	uint16_t TDC_STAT_CARD_B_FLGS = device("GPIO_SW_LEDS")->readRegister("TDC_STAT_CARD_B_FLGS");

	// See email from G. Wittwer, March 2016
	// In external clock mode, GMC_SEL needs to be configured; typically for a slave MuTanT module connected in standard VHDCI GMC_SEL<1..0> = 2.
	// This could be the default value.
	utl::BitFieldHelper< uint16_t >::setField(TDC_STAT_CARD_B_FLGS, 8, 2, config.externalGmcClockSelection);

	NODE_DEBUG() << "Setting register 'TDC_STAT_CARD_B_FLGS' to " << std::hex << std::showbase << TDC_STAT_CARD_B_FLGS << std::dec;
	device("GPIO_SW_LEDS")->writeRegister("TDC_STAT_CARD_B_FLGS", TDC_STAT_CARD_B_FLGS);
}
//_________________________________________________________________________________________________
/**
 * Configures registers of MuTanT B board.
 */
void MutantNode::configureRegistersMutantB(CCfg::CConfig /* cfg */)
{
	NODE_DEBUG() << "Configuring MuTanT B registers";
	try
	{
		// Select clock source
		selectClockSource(config.isSlave, config.isSlavePortA);
		// Set register GLOBAL_CONF
		setGlobalConfigurationRegister();
		// Enable test mode
		selectTestMode();
		// Enable readout mode
		selectReadoutMode();
		// Set timestamp and event counters
		setTimestampAndEventCounters();
		// Initialize scalers
		configureScalers();
		// Set source ID
		setSourceId(config.sourceId);
		// Select frame type
		setFrameType(config.frameType);
		// Set shelf mask
		setShelfMask();
		// Set AMC_TAG_ONLY
		setAmcTagOnly();
		// Set user register
		setUserRegister();
		// Configure Level 0
		configureLevel0Trigger();
		// Configure Level 1
		configureLevel1Trigger();
		// Configure logical inspection signals
		configureLogicalInspectionSignals();
		// Check status
		checkStatus();
		checkAmcStatus();
	}
	catch (const CCfg::CConfig::NotFound & e)
	{
		throw CmdException(e.what());
	}
}
//_________________________________________________________________________________________________
/**
 * Register Beast boards in MuTanT.
 *
 * We assume here that any board declared in AMC_AG_ONLY should be declared in CLOCK_TREE_CONFIG1.
 */
void MutantNode::registerBeastBoards(CCfg::CConfig /* cfg */)
{
	if (config.amcTagOnly <= 0)
		return;

	NODE_DEBUG() << "Updating CLOCK_TREE_CONFIG1 to take into account BEAST or other non-CoBo boards" << std::hex << std::showbase << config.amcTagOnly << std::dec;
	try
	{
		// Read register CLOCK_TREE_CONFIG1
		uint16_t amc_clk1_out_en  = device("GPIO_CLK_TREE")->readField("CLOCK_TREE_CONFIG1", "amc_clk1_out_en");
		uint16_t amc_clk2_out_en = device("GPIO_CLK_TREE")->readField("CLOCK_TREE_CONFIG1", "amc_clk2_out_en");
		// Edit bit pattern to include BEAST boards declared in AMC_TAG_ONLY
		for (size_t slot=0; slot<12; slot++)
		{
			bool ato = utl::BitFieldHelper< uint16_t >::getField(config.amcTagOnly, slot, 1u);
			if (ato)
			{
				utl::BitFieldHelper< uint16_t >::setField(amc_clk1_out_en, slot, 1u, ato);
				utl::BitFieldHelper< uint16_t >::setField(amc_clk2_out_en, slot, 1u, ato);
			}
		}
		// Modify register
		setClockShelfMasks(amc_clk1_out_en, amc_clk2_out_en);
	}
	catch (const mdaq::utl::CmdException &)
	{
		NODE_WARN() << "Error registering BEAST boards in CLOCK_TREE_CONFIG1.";
		throw;
	}
}
//_________________________________________________________________________________________________
/**
 * Set 8 bit  SPI digital potentiometer tuning the frequency of the VCXO chip.
 * @param vcxoControl 8-bit value
 *
 * Email from G. Wittwer on Jan 19, 2017:
 * I used also a VCXO on MUTANT which has a frequency of 200 MHz, a frequency  stability vs temperature of 10 ppm
 * and an available pulling range of +/- 100 ppm (that is not a tolerance) in order to tune this frequency easily
 * if tomorrow MUTANT must be synchronized with other clock module at GANIL where F0 is never exactly 100 000 000 Hz.
 *
 * So, the VC (Voltage Control pin) of the VCXO used on MUTANT is driven by an 8 bit  SPI digital potentiometer for tuning.
 * Normally this one is initialized internally, at power up, to what they called MID-SCALE (0x80) but anyway,
 * I initialize it again (by default) in my standalone programs to 0x80 and I obtain most of the time a value of 99 998 300 Hz
 * that is only 1.7KHz from 100 000 000 Hz. (I use a 10 digit frequency meter able to display the Hertz for all my measures)
 *
 * You can program a value, just a bit higher, in order to adjust it to few tenth of Hertz but I never need to do that.
 * In fact it seems at power up that when MUTANT is booting LINUX, the digi-pot is set to 0 but then this one is never initialized.
 * So in this case the output frequency, corresponds to the lowest frequency value, that you measured: 99.9879MHz
 * typically 12100 Hz below 100 000 000 Hz, of course out of the +/- 5KHz range required by CoBo PLL.
 */
void MutantNode::tuneVcxoFrequency(const uint8_t vcxoControl)
{
	NODE_DEBUG() << "Tuning VCXO frequency with control set to " <<  std::hex << std::showbase << (uint16_t) vcxoControl << std::dec;
	device("CLK200_SPI")->writeRegister("VCXO_CONTROL", vcxoControl);
}
//_________________________________________________________________________________________________
/**
 * Selects the clock source used by the trigger/timestamper (FPGA board B)
 * and manages the WSCA signal as an output if the current MUTANT is master or as an input if it is slave.
 *
 * Sets the GMC_WSCA_SEL register of device CARD_B_EMC
 *
 * @see MuTanT Slow Control Registers (v1.5)
 */
void MutantNode::selectClockSource(const bool & isSlave, const bool & slavePortA)
{
	// Set GMC source in register GMC_WSCA_SEL of device CARD_B_EMC
	// New table in slow control register mapping v1.5
	// GMC  = 2: Directly from MUTANT_A VCXO
	//      = 0: From MUTANT_A but phase aligned (same path for Master or Slave)
	//      = 1: From the optical link of MUTANT_A
	//      = 3: MUTANT_B 100 MHz XO (for MUTANT_B standalone tests)
	// Value = 0 = From MUTANT_A but phase aligned for master or slave
	const uint8_t gmc = 0;
	device("CARD_B_EMC")->writeField("GMC_WSCA_SEL", "GMC", gmc);

	// Set WSCA source in GMC_WSCA_SEL register of device CARD_B_EMC
	// WSCA<1..0> = 0: LOCAL if (master)
	//            = 1: REMOTE from VHDCI port A (with slave MUTANT)
	//            = 2: REMOTE from VHDCI port B (with slave MUTANT)
	//            = 3: Not Used
	uint8_t wsca = 0;
	if (isSlave)
	{
		wsca = slavePortA?UINT8_C(1):UINT8_C(2);
	}
	device("CARD_B_EMC")->writeField("GMC_WSCA_SEL", "WSCA", wsca);

	// Set WSCA output path at the MASTER level (VHDCI connector – only used for slave)
	// WSCA<2> = 0: Delayed WSCA output (managed)
	//         = 1: Prompt WSCA (no delay management)
	// Always set to 0 in slave, set to 1 in master if there is at least 1 slave
	bool promptWsca = false;
	if (mutantModuleCount_ >= 2 and not config.isSlave)
		promptWsca = true;
	device("CARD_B_EMC")->writeField("GMC_WSCA_SEL", "PromptWSCA", promptWsca);

	NODE_DEBUG() << "Selecting clock source GMC=" << (uint16_t) gmc << " WSCA=" << (uint16_t) wsca << " PromptWSCA=" << (short) promptWsca;
}
//_________________________________________________________________________________________________
/**
 * Sets the register GLOBAL_CONF of device CARD_B_EMC
 *
 * @see MuTanT Slow Control Registers (v1.4)
 */
void MutantNode::setGlobalConfigurationRegister()
{
	// Build register value
	uint16_t regValue = 0;

	// Bit 0: M/S = 1: this MUTANT module is master
	//            = 0: this MUTANT module is slave
	utl::BitFieldHelper< uint16_t >::setField(regValue, 0, 1, config.isSlave?0:1);

	// Bit 1: SNUM: if the current MUTANT is slave, this bit gives the slave number (0 or 1).
	// It also enables the connection to the right VHDCI port (A or B).
	// This bit has no meaning in case of master MUTANT
	utl::BitFieldHelper< uint16_t >::setField(regValue, 1, 1, config.slaveId);

	// Bit 2: VPA_ON = 1 VHDCI port A is enabled. A Master or slave MUTANT is connected via this port
	//        VPA_ON = 0 VHDCI port A is disabled
	utl::BitFieldHelper< uint16_t >::setField(regValue, 2, 1, config.vhdciPortAEnabled);

	// Bit 3: VPB_ON = 1 VHDCI port B is enabled. A Master or slave MUTANT is connected via this port
	//        VPB_ON = 0 VHDCI port B is disabled
	utl::BitFieldHelper< uint16_t >::setField(regValue, 3, 1, config.vhdciPortBEnabled);

	// Bit 4: L0EN = 1 -> trigger Level 0 (external) is enabled
	utl::BitFieldHelper< uint16_t >::setField(regValue, 4, 1, config.isSlave?0:config.level0Enabled);

	// Bit 5: L1EN = 1 -> trigger Level 1 (multiplicity) is enabled
	utl::BitFieldHelper< uint16_t >::setField(regValue, 5, 1, config.isSlave?0:config.level1Enabled);

	// Bit 6: L2EN = 1 -> trigger Level 2 (hit pattern) is enabled
	utl::BitFieldHelper< uint16_t >::setField(regValue, 6, 1, config.isSlave?0:config.level2Enabled);

	// Bit 7: L2W = 1 tells the trigger/timestamper state machine that a mask pattern will be provided to CoBos in case of L2 OK
	utl::BitFieldHelper< uint16_t >::setField(regValue, 7, 1, config.level2Enabled ? config.hitPatternModified : 0);

	// Bit 8: L0/L1 = 0: External trigger (L 0) starts a coincidence window waiting for a valid (L1) event (TPC multiplicity)
	//              = 1: When Level 1 trigger occurs, a coincidence window is opened waiting for a L0 trigger validation
	// This bit has no effect unless both L0 and L1 are enabled.
	utl::BitFieldHelper< uint16_t >::setField(regValue, 8, 1, config.l0l1);

	// Doc 1.4: Bits 9-10: L1SEL<1..0> = 00: Logic “OR” of the 2 multiplicity thresholds for a valid Level 1 OK
	//                                 = 01: Only the low multiplicity threshold is taken into account for a valid Level 1 OK
	//                                 = 10: Only the high multiplicity threshold is taken into account for a valid Level 1 OK
	//                                 = 11: Logic “AND” of the 2 multiplicity thresholds for a valid Level 1 OK
	// Doc 1.5: Bits 9-10: L1SEL<1..0> = 00: Only the low multiplicity threshold is taken into account for a valid Level 1 OK
	//                                 = 01: Only the high multiplicity threshold is taken into account for a valid Level 1 OK
	//                                 = 10: Logic “OR” of the 2 multiplicity thresholds for a valid Level 1 OK
	//                                 = 11: Logic “AND” of the 2 multiplicity thresholds for a valid Level 1 OK
	utl::BitFieldHelper< uint16_t >::setField(regValue, 9, 2, config.level1Selection);

	// Bit 11: TPD = 0: Normal mode
	//             = 1: 2 proton decay mode (double set of registers is used)
	utl::BitFieldHelper< uint16_t >::setField(regValue, 11, 1, config.tpdModeEnabled);

	// Bit 13: HER = 0: Normal mode
	//             = 1: Half-Event Readout (in 2p mode, half-event sent event after timeout)
	utl::BitFieldHelper< uint16_t >::setField(regValue, 12, 1, config.halfEventReadoutEnabled);

	// Bit 14: C/S = 0: The CENTRUM/SCALERS interface is disabled
	//             = 1: The CENTRUM/SCALERS interface is enabled
	// Bit 15: CENT = 0: The scaler channels are operational if C/S bit is one
	//              = 1: The coupling with CENTRUM is enabled if C/S bit is one
	switch (config.csMode)
	{
	case Scalers:
		utl::BitFieldHelper< uint16_t >::setField(regValue, 14, 1, 1);
		utl::BitFieldHelper< uint16_t >::setField(regValue, 15, 1, 0);
		break;
	case CENTRUM:
		utl::BitFieldHelper< uint16_t >::setField(regValue, 14, 1, 1);
		utl::BitFieldHelper< uint16_t >::setField(regValue, 15, 1, 1);
		break;
	default:
		utl::BitFieldHelper< uint16_t >::setField(regValue, 14, 1, 0);
		utl::BitFieldHelper< uint16_t >::setField(regValue, 15, 1, 0);
	}

	// Set register
	NODE_DEBUG() << "Setting global configuration register to " << std::hex << std::showbase << regValue << std::dec;
	device("CARD_B_EMC")->writeRegister("GLOBAL_CONF", regValue);
}
//_________________________________________________________________________________________________
/**
 * BR 13 may 2015
 * Enables or disables readout mode by setting the bit READOUT of register TRIGGER_CTRL of device CARD_B_EMC.
 */
void MutantNode::selectReadoutMode()
{
	if (config.readoutModeEnabled)
	{
		// Set field
		device("CARD_B_EMC")->writeField("TRIGGER_CTRL", "read_mode", 0);
	}
	else
	{
		NODE_WARN() << "Disabling readout mode of Mutant[" << moduleId_ << ']';
		// Set field
		device("CARD_B_EMC")->writeField("TRIGGER_CTRL", "read_mode", 1);
	}
}
//_________________________________________________________________________________________________
/**
 * Enables or disables test mode by setting the bit TEST of register TRIGGER_CTRL of device CARD_B_EMC.
 */
void MutantNode::selectTestMode()
{
	if (config.testModeEnabled)
	{
		NODE_WARN() << "Enabling test mode of Mutant[" << moduleId_ << ']';
	}
	// Set field
	device("CARD_B_EMC")->writeField("TRIGGER_CTRL", "test", config.testModeEnabled);
}
//_________________________________________________________________________________________________
/**
 * Sets timestamp and event counters
 *
 * Register TRIGGER_CMD: CL2C CL1BC CL1AC CLØC EEC LEC EUC LUC
 */
void MutantNode::setTimestampAndEventCounters()
{
	NODE_DEBUG() << "Initializing timestamp and event counters";
	// Set EUC and EEC bits to 0 to disable Universal Clock counter and Event counter
	// EUC = 1: Enable Universal Clock counter
	// EEC = 1: Enable Event Counter
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "euc", 0);
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "eec", 0);

	// Clear event counters
	if (config.clearCountersEnabled)
	{
		clearEventCounters();
	}

	// Set SET_TIME_xx and SET_EVT_NUM_xxx registers to user-defined values

	if (config.preLoadTimestamp)
	{
		NODE_DEBUG() << "Pre-loading timestamp " << config.timestamp_10ns << " x 10 ns";
		// SET_TIME_xxx: MUTANT 48 bit timestamp counter (in units of 10 ns) in three 16-bit chunks
		device("CARD_B_EMC")->writeRegister("SET_TIME_LOW",  utl::BitFieldHelper< uint64_t >::getField(config.timestamp_10ns,  0u, 16u));
		device("CARD_B_EMC")->writeRegister("SET_TIME_MID",  utl::BitFieldHelper< uint64_t >::getField(config.timestamp_10ns, 16u, 16u));
		device("CARD_B_EMC")->writeRegister("SET_TIME_HIGH", utl::BitFieldHelper< uint64_t >::getField(config.timestamp_10ns, 32u, 16u));

		// Set LUC bit to 1
		// LUC = 1: Load Universal Clock counter with the value preloaded in the SET_TIME_xxx registers
		device("CARD_B_EMC")->writeField("TRIGGER_CMD", "luc", 1);
		// Set LUC bit back to 0
		device("CARD_B_EMC")->writeField("TRIGGER_CMD", "luc", 0);

		// Set L2PTC to 1 and back to 0
		// L2PTC= 1: Load "2 proton" 32 bit time counter with SET_TIME_MID & SET_TIME_LOW register values
		try
		{
			device("CARD_B_EMC")->writeField("TRIGGER_CMD", "l2ptc", 1);
			device("CARD_B_EMC")->writeField("TRIGGER_CMD", "l2ptc", 0);
		}
		catch (const mdaq::utl::CmdException & e)
		{
			NODE_WARN() << e.reason;
			NODE_WARN() << "You should update your MuTanT description file";
		}
	}

	if (config.preLoadEventNumber)
	{
		NODE_DEBUG() << "Pre-loading event number " << config.eventNumber;
		// SET_EVT_NUM_xxx: MUTANT 32 bit event counter in two 16-bit chunks
		device("CARD_B_EMC")->writeRegister("SET_EVT_NUM_LOW",  utl::BitFieldHelper< uint64_t >::getField(config.eventNumber,  0u, 16u));
		device("CARD_B_EMC")->writeRegister("SET_EVT_NUM_HIGH", utl::BitFieldHelper< uint64_t >::getField(config.eventNumber, 16u, 16u));

		// Set LEC bit to 1
		// LEC = 1: Load Event Counter with the value preloaded in the SET_EVT_NUM_xxx registers
		device("CARD_B_EMC")->writeField("TRIGGER_CMD", "lec", 1);
		// Set LEC bit back to 0
		device("CARD_B_EMC")->writeField("TRIGGER_CMD", "lec", 0);
	}

	// Set EUC and EEC bits to 1 to enable Universal Clock counter and Event counter
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "euc", 1);
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "eec", 1);
}
//_________________________________________________________________________________________________
/**
 * Enables or disables counting for a given scaler.
 *
 * Sets the ENS? bit of register SCALER_CTRL of device CARD_B_EMC.
 *
 * @param scalerId Scaler number.
 * @param enable Whether to enable or disable the scaler.
 */
void MutantNode::enableScaler(const size_t & scalerId, const bool enable)
{
	// ENS1= 1: Enable counting for scaler 1
	// ENS2= 1: Enable counting for scaler 2
	// ENS3= 1: Enable counting for scaler 3
	// ENS4= 1: Enable counting for scaler 4
	// ENS5= 1: Enable counting for scaler 5

	const std::string fieldName = std::string("ens") + lexical_cast< std::string >(scalerId);
	device("CARD_B_EMC")->writeField("SCALER_CTRL", fieldName, enable?UINT8_C(1):0);
}
//_________________________________________________________________________________________________
/**
 * Enables or disables counting for scalers.
 *
 * Sets the ENS? bits of register SCALER_CTRL of device CARD_B_EMC.
 */
void MutantNode::enableScalers(bool enable)
{
	// ENS1= 1: Enable counting for scaler 1
	// ENS2= 1: Enable counting for scaler 2
	// ENS3= 1: Enable counting for scaler 3
	// ENS4= 1: Enable counting for scaler 4
	// ENS5= 1: Enable counting for scaler 5

	for (size_t i=1; i <=5; ++i)
	{
		enableScaler(i, enable);
	}
}
//_________________________________________________________________________________________________
/**
 * Clears event counters.
 *
 * Sets the CL?C bits of register TRIGGER_CMD of device CARD_B_EMC to 1 and back to 0.
 */
void MutantNode::clearEventCounters()
{
	// Set CL*C bits to 1 and back to 0 to clear the accepted event counters
	// CL0C= 1: Clear Level 0 accepted event counter
	// CL1AC= 1: Clear Level 1A accepted event counter
	// CL1BC= 1: Clear Level 1B accepted event counter
	// CL2C= 1: Clear Level 2 accepted event counter
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "cl0c", 1);
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "cl1ac", 1);
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "cl1bc", 1);
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "cl2c", 1);

	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "cl0c", 0);
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "cl1ac", 0);
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "cl1bc", 0);
	device("CARD_B_EMC")->writeField("TRIGGER_CMD", "cl2c", 0);
}
//_________________________________________________________________________________________________
/**
 * Clears scalers.
 *
 * Sets the CS? bits of register SCALER_CTRL of device CARD_B_EMC to 1 and back to 0.
 */
void MutantNode::clearScalers()
{
	// CS1= 1: Clear scaler 1
	// CS2= 1: Clear scaler 2
	// CS3= 1: Clear scaler 3
	// CS4= 1: Clear scaler 4
	// CS5= 1: Clear scaler 5

	for (size_t i=1; i <=5; ++i)
	{
		const std::string fieldName = std::string("cs") + lexical_cast< std::string >(i);
		device("CARD_B_EMC")->writeField("SCALER_CTRL", fieldName, 1);
		device("CARD_B_EMC")->writeField("SCALER_CTRL", fieldName, 0);
	}
}
//_________________________________________________________________________________________________
/**
 * Initializes scalers. Sets the register SCALER_CTRL of device CARD_B_EMC.
 *   - Disables all scalers
 *   - Clears scaler counters if required by configuration.
 *   - Enables those of the scalers selected in the configuration.
 */
void MutantNode::configureScalers()
{
	NODE_DEBUG() << "Initializing scalers";
	// Disable counting for all scalers
	enableScalers(false);

	// Clear scalers
	if (config.clearScalersEnabled)
	{
		clearScalers();
	}

	// Re-enable scalers if required
	for (size_t i=1; i <=5; ++i)
	{
		if (config.enableScalers[i])
			enableScaler(i, true);
	}
}
//_________________________________________________________________________________________________
/**
 * Sets the register SOURCE_SYSTEM of device CARD_B_EMC.
 *
 * The first byte of this register will be included in the readout of the data frames coming from the FIFO
 *  in order to identify the detection system of which MUTANT is part.
 */
void MutantNode::setSourceId(const uint8_t sourceId)
{
	NODE_DEBUG() << "Setting source ID to " << (uint16_t) sourceId;
	device("CARD_B_EMC")->writeField("SOURCE_SYSTEM", "source", sourceId);
}
//_________________________________________________________________________________________________
/**
 * Sets the register FRAME_TYPE of device CARD_B_EMC.
 *
 * FT<0..1> = 00: A basic frame is built on each accepted event (see MUTANT-DATA readout manual)
 *          = 01: Basic frame version + trigger counters are used to built a frame on each accepted event
 *          = 10: This frame is based on the previous one plus the ten words of the five scalers
 *          = 11: Reserved
 */
void MutantNode::setFrameType(const uint8_t type)
{
	NODE_DEBUG() << "Selecting frame type no. " << (uint16_t) type;
	device("CARD_B_EMC")->writeField("FIFO_FRAME_TYPE", "ft", type);
}
//_________________________________________________________________________________________________
/**
 * Save value of CLOCK_TREE_CONFIG1 after mono-crate alignment,
 *  as it will be overridden by multi-crate alignment
 *  and will be needed later on to set shelf mask.
 */
void MutantNode::storeClockShelfMasks()
{
	// Read register CLOCK_TREE_CONFIG1
	fMaskGmc_  = device("GPIO_CLK_TREE")->readField("CLOCK_TREE_CONFIG1", "amc_clk1_out_en");
	fMaskWsca_ = device("GPIO_CLK_TREE")->readField("CLOCK_TREE_CONFIG1", "amc_clk2_out_en");
}
//_________________________________________________________________________________________________
/**
 * Set enabled GMC and WSCA outputs (register CLOCK_TREE_CONFIG1)
 */
void MutantNode::setClockShelfMasks(const uint16_t & gmcMask, const uint16_t & wscaMask)
{
	// Read register CLOCK_TREE_CONFIG1
	device("GPIO_CLK_TREE")->writeField("CLOCK_TREE_CONFIG1", "amc_clk1_out_en", gmcMask);
	device("GPIO_CLK_TREE")->writeField("CLOCK_TREE_CONFIG1", "amc_clk2_out_en", wscaMask);
}
//_________________________________________________________________________________________________
/**
 * Restore enabled GMC and WSCA outputs (register CLOCK_TREE_CONFIG1)
 */
void MutantNode::restoreClockShelfMasks()
{
	setClockShelfMasks(fMaskGmc_, fMaskWsca_);
}
//_________________________________________________________________________________________________
/**
 * Sets the register SHELF_MASK of device CARD_B_EMC.
 *
 * Bits 0 to 11 indicate whether the CoBo boards in AMC slots 1 to 12 are enabled (bit set) or disabled.
 * The mask is set based on
 *   - the contents of register CLOCK_TREE_CONFIG1 of MuTanT A device GPIO_CLOCK_TREE.
 *   - the "isActive" CoBo user parameters.
 *
 * On April 30th, 2015:
 * - the register MASTER_MASK was renamed into SHELF_MASK.
 * - the unused registers SLAVE_1_MASK and SLAVE_2_MASK have been renamed into USER_REG1 and USER_REG2.
 *
 */
void MutantNode::setShelfMask()
{
	uint16_t master_mask = UINT16_C(0);

	// Loop over AMC slots
	for (size_t amc=0; amc<12u; ++amc)
	{
		// Check both clocks are OK
		if (not (utl::BitFieldHelper< uint16_t >::getBit(fMaskGmc_, amc) and utl::BitFieldHelper< uint16_t >::getBit(fMaskWsca_, amc)))
			continue;

		// Enable slot
		utl::BitFieldHelper< uint16_t >::setBit(master_mask, amc, true);
	}
	// Set master mask
	NODE_DEBUG() << "Setting shelf mask to " << std::hex << std::showbase << master_mask << std::dec;
	try
	{
		device("CARD_B_EMC")->writeRegister("SHELF_MASK", master_mask);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		NODE_WARN() << "Error setting register SHELF_MASK. Trying its old name MASTER_MASK instead.";
		device("CARD_B_EMC")->writeRegister("MASTER_MASK", master_mask);
	}
}
//_________________________________________________________________________________________________
/**
 * Set ATO bits e.g. define AMC slots with board (typically BEAST) that should receive only synchronization information.
 */
void MutantNode::setAmcTagOnly()
{
	NODE_DEBUG() << "Setting AMC_TAG_ONLY register to " << std::hex << std::showbase << config.amcTagOnly << std::dec;
	try
	{
		device("CARD_B_EMC")->writeRegister("AMC_TAG_ONLY", config.amcTagOnly);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		NODE_WARN() << "Error setting register AMC_TAG_ONLY. Trying its old name USER_REG1 instead.";
		device("CARD_B_EMC")->writeRegister("USER_REG1", config.amcTagOnly);
	}
}
//_________________________________________________________________________________________________
/**
 * Set general purpose register not attached to hardware where any information can be stored.
 */
void MutantNode::setUserRegister()
{
	NODE_DEBUG() << "Setting user register to " << std::hex << std::showbase << config.userReg << std::dec;
	try
	{
		device("CARD_B_EMC")->writeRegister("USER_REG", config.userReg);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		NODE_WARN() << "Error setting register USER_REG. Trying its old name USER_REG2 instead.";
		device("CARD_B_EMC")->writeRegister("USER_REG2", config.userReg);
	}
}
//_________________________________________________________________________________________________
/**
 * Sets the registers related to Level 0 trigger.
 */
void MutantNode::configureLevel0Trigger()
{
	NODE_DEBUG() << "Configuring Level 0 trigger";

	// Set register L0_GATE_DELAY (delay of the Level 0 gate, in units of 10 ns, on 16 bits)
	device("CARD_B_EMC")->writeRegister("L0_GATE_DELAY", config.level0GateDelay_10ns);

	// Set register L0_GATE_DURATION (duration of the Level 0 gate, in units of 10 ns, on 16 bits)
	device("CARD_B_EMC")->writeRegister("L0_GATE_DURATION", config.level0GateDuration_10ns);
}
//_________________________________________________________________________________________________
/**
 * Sets the registers related to the Level 1 trigger.
 */
void MutantNode::configureLevel1Trigger()
{
	NODE_DEBUG() << "Configuring Level 1 trigger";
	RemoteRegisterProcessor & regProc = hwNode().regProc("CARD_B_EMC");

	regProc
		// Set low multiplicity threshold of Level 1A trigger: register L1_A_LOW_MULT_THRESHOLD
		.write("L1_A_LOW_MULT_THRESHOLD", config.level1aLowThreshold)

		// Set high multiplicity threshold of Level 1A trigger: register L1_A_HIGH_MULT_THRESHOLD
		.write("L1_A_HIGH_MULT_THRESHOLD", config.level1aHighThreshold)

		// Set low multiplicity threshold of Level 1B trigger: register L1_B_LOW_MULT_THRESHOLD
		.write("L1_B_LOW_MULT_THRESHOLD", config.level1bLowThreshold)

		// Set high multiplicity threshold of Level 1B trigger: register L1_B_HIGH_MULT_THRESHOLD
		.write("L1_B_HIGH_MULT_THRESHOLD", config.level1bHighThreshold)

		// Set number of buckets for Level 1A low trigger: register L1_A_LOW_TH_N_BUCKETS
		.write("L1_A_LOW_TH_N_BUCKETS", config.level1aLowNumBuckets)

		// Set number of buckets for Level 1A high trigger: register L1_A_HIGH_TH_N_BUCKETS
		.write("L1_A_HIGH_TH_N_BUCKETS", config.level1aHighNumBuckets)

		// Set number of buckets for Level 1B low trigger: register L1_B_LOW_TH_N_BUCKETS
		.write("L1_B_LOW_TH_N_BUCKETS", config.level1bLowNumBuckets)

		// Set number of buckets for Level 1B high trigger: register L1_B_HIGH_TH_N_BUCKETS
		.write("L1_B_HIGH_TH_N_BUCKETS", config.level1bHighNumBuckets)

		// Set timeout value in case of no L1B trigger: register L1_B_TIMEOUT
		// Time base in v1.4.0: LSB =10 ns, min value = 10 ns, max value = 655.35 μs
		// Time base in v1.4.1: LSB =1 μs, min value = 0, max value = 65.535 ms
		// Time base in v1.5.0: LSB = 2,56 μs, min value = 1, max value ≈ 167 ms
		.write("L1_B_TIMEOUT", config.level1bTimeout_2560ns)

		// Set time added to a valid L1 trigger in order to delay the STOP sampling: register L1_POST_TRIG_DELAY
		.write("L1_POST_TRIG_DELAY", config.level1PostTrigDelay_10ns);

	if (config.level1bTimeout_2560ns < 1u)
	{
		NODE_WARN() << "Setting L1B timeout to 0 us!";
	}
}
//_________________________________________________________________________________________________
/**
 * Selects the internal signals sent to the front panel LEMO connectors "LOG INSP1" and "LOG INSP2"
 *  and the internal signals that will be sent on the rear panel LEMO connectors "TEST1" and "TEST2".
 */
void MutantNode::configureLogicalInspectionSignals()
{
	uint16_t regValue = 0;

	// Select LOG INSP1 signal
	NODE_DEBUG() << "Selecting signal " << (uint16_t) config.frontPanelInsp1 << " for front panel LEMO 'LOG INSP1'";
	utl::BitFieldHelper< uint16_t >::setField(regValue, 0, 8u, config.frontPanelInsp1);

	// Select LOG INSP2 signal
	NODE_DEBUG() << "Selecting signal " << (uint16_t) config.frontPanelInsp2 << " for front panel LEMO 'LOG INSP2'";
	utl::BitFieldHelper< uint16_t >::setField(regValue, 8u, 8u, config.frontPanelInsp2);

	// Set register FP_LOG_IREG
	device("CARD_B_EMC")->writeRegister("FP_LOG_IREG", regValue);

	regValue = 0;

	// Select TEST1 signal
	NODE_DEBUG() << "Selecting signal " << (uint16_t) config.rearPanelInsp1 << " for rear panel LEMO 'TEST1'";
	utl::BitFieldHelper< uint16_t >::setField(regValue, 0, 8u, config.rearPanelInsp1);

	// Select TEST2 signal
	NODE_DEBUG() << "Selecting signal " << (uint16_t) config.rearPanelInsp2 << " for rear panel LEMO 'TEST2'";
	utl::BitFieldHelper< uint16_t >::setField(regValue, 8u, 8u, config.rearPanelInsp2);

	// Set register RP_LOG_IREG
	device("CARD_B_EMC")->writeRegister("RP_LOG_IREG", regValue);
}
//_________________________________________________________________________________________________
/**
 * Checks contents of the MUTANT_STATUS register.
 * Check master shelf and slave crates are ready.
 */
void MutantNode::checkStatus()
{
	const uint16_t status = device("CARD_B_EMC")->readRegister("MUTANT_STATUS");
	NODE_DEBUG() << "Checking MuTanT status register: " << std::hex << std::showbase << status << std::dec << std::noshowbase;

	// For bit MRDY to be set,  bits 'mutantMode' and 'initDone' must be set in register 'mutantConfig' of CoBo.
	const bool masterReady = device("CARD_B_EMC")->readField("MUTANT_STATUS", "mrdy");
	if (not masterReady)
	{
		throw CmdException("Status register indicates master shelf is not ready!");
	}
	/* BR. 2014/10/08 tester si esclave enabled avant de tester ces bits de status
	const bool s0Ready = device("CARD_B_EMC")->readField("MUTANT_STATUS", "s0rdy");
	if (not s0Ready)
	{
		throw CmdException("Status register indicates not all CoBo boards in slave 0 crate are ready!");
	}

	const bool s1Ready = device("CARD_B_EMC")->readField("MUTANT_STATUS", "s1rdy");
	if (not s1Ready)
	{
		throw CmdException("Status register indicates not all CoBo boards in slave 1 crate are ready!");
	}
	*/
}
//_________________________________________________________________________________________________
/**
 * Checks status of AMC slots
 * Reads registers AMCx_ySTATUS
 *
 * @see Communication protocol with the CoBo boards (V1.1)
 */
void MutantNode::checkAmcStatus()
{
	// Loop over AMC
	for (size_t amc=1u; amc <= 12u; ++amc)
	{
		// Read current state code of AMC on 4 bits
		const size_t regId = (amc - 1u)/4u;
		std::string regName = boost::str(boost::format("AMC%1%_%2%_STATUS") % (4*regId+4) % (4*regId+1));
		std::string fieldName = boost::str(boost::format("amc%1%_status_code") % amc);
		const uint16_t code = device("CARD_B_EMC")->readField(regName, fieldName);
		NODE_DEBUG() << "Status of AMC " << amc << ": " << code;
	}
}
//_________________________________________________________________________________________________
/**
 * Set module is running by setting the RUN bit of the TRIGGER_CTRL register of MuTanT B device CARB_B_EMC.
 *
 * BR. 2014/10/09
 */
void MutantNode::setModuleRunning()
{
	if (isTriggerEnabled() or isSlaveModule())
	{
		NODE_DEBUG() << "Setting MuTanT 'run' bit!";
		device("CARD_B_EMC")->writeField("TRIGGER_CTRL", "run", 1);
	}
	else
	{
		NODE_WARN() << "No trigger mode was enabled, MuTanT 'run' bit WILL NOT be set!";
	}
}
//_________________________________________________________________________________________________
/**
 * Performs pre-run actions.
 */
void MutantNode::runInit()
{
	// Clear scaler counters
	if (config.clearScalersEnabled)
	{
		clearScalers();
	}
	// Clear event counters
	if (config.clearCountersEnabled)
	{
		clearEventCounters();
	}
}
//_________________________________________________________________________________________________
void MutantNode::impl_start()
{
	if (not isSlaveModule())
	{
		if (isTriggerEnabled())
		{
			if (isReadoutEnabled())
			{
				if (hwNode().dataSenderCtrlPrx)
				{
					hwNode().dataSenderCtrlPrx->daqStart();
				}
			}
			else
			{
				NODE_INFO() << "MuTanT readout is disabled";
			}

			// Set Module running (BR. 2014/10/09)
			setModuleRunning();
		}
		else
		{
			NODE_WARN() << "No trigger mode was enabled, MuTanT 'run' bit WILL NOT be set!";
		}
	}
	else
	{
		// Set Module running (BR. 2014/10/09)
		setModuleRunning();
	}
}
//_________________________________________________________________________________________________
void MutantNode::impl_stop()
{
	if (isReadoutEnabled() and not isSlaveModule())
	{
		if (hwNode().dataSenderCtrlPrx)
		{
			hwNode().dataSenderCtrlPrx->daqStop();
		}
    	}

	// Set Module stopped (BR. 2014/10/09)
	ensureModuleNotRunning();
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
