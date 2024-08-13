/**
 * @file SystemManager.cpp
 * @date Mar 18, 2014
 * @author sizun
 *
 * @note SVN tag: $Id: SystemManager.cpp 1739 2017-11-23 16:25:34Z psizun $
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

#include "get/rc/SystemManager.h"
#include "get/mutant/AlignmentServer.h"
#include "get/rc/AmcSlotManager.h"
#include "get/rc/BackEnd.h"
#include "get/rc/NodeId.h"
#include <string>

using ::mdaq::utl::CmdException;

namespace get
{
namespace rc
{

/**
 * Constructor creates the 2 list of MuTanT and CoBo boards.
 * @param backend Reference to backend callbacks.
 */
SystemManager::SystemManager(BackEnd& backend) : Node::ListList(backend)
{
	addNewList<CoBoNode>();
	addNewList<MutantNode>();
}

/**
 * De-allocates all resources that are not automatically disposed of.
 */
SystemManager::~SystemManager()
{
}

/**
 * Destroys all nodes associated to electronics boards (CoBo, Mutant...)
 */
void SystemManager::removeNodes()
{
	Node::ListList::removeNodes();
	amcSlots.clear();
}

CoBoNode & SystemManager::getCoBoNode(const std::string & name)
{
	for (Node::List::Iterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		if (cobo->hwNode().hwServerCtrlPrx->name() == name)
			return *cobo;
	}
	throw mdaq::utl::CmdException(std::string("Error accessing CoBoNode instance named '") + name + "'");
}

/**
 * Returns the Master MuTanT module.
 * @return
 */
MutantNode & SystemManager::mutantMaster()
{
	for (Node::List::Iterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		if (not mutant->isSlaveModule())
			return *mutant;
	}
	throw CmdException("Could not find any master MuTanT module!");
}

/**
 * Returns the Master MuTanT module.
 * @return
 */
bool SystemManager::hasMutantMaster()
{
	try
	{
		mutantMaster();
	}
	catch (const CmdException &)
	{
		return false;
	}
	return true;
}

MutantSeq SystemManager::mutantSlaves()
{
	MutantSeq slaves;
	for (Node::List::Iterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		if (mutant->isSlaveModule())
			slaves.push_back(&(*mutant));
	}
	return slaves;
}

void SystemManager::describe()
{
	LOG_INFO() << "Connected to "
		       << getList<CoBoNode>().size() << " CoBo board(s) and "
		       << getList<MutantNode>().size() << " MuTanT module(s)";
	// Make all CoBo boards send <em>not configured</em> status to MuTanT modules.

	setCoBoInitDone(false);
	// BR. 2014/10/08 doesn't work, must be done at the end of the configure phase
	// Tell CoBo boards whether they will be MuTanT driven
	//	setCoBoMutantMode(mutantManagers.size() > 0);
	setCoBoMutantMode(false);
}

/**
 * Sets @em initDone bit of all CoBo boards to 0 or 1 so that they send <em>not configured</em> or <em>ready</em> status to Mutant.
 */
void SystemManager::setCoBoInitDone(bool done)
{
	if (done)
	{
		LOG_INFO() << "Setting 'initDone' bit of CoBo board(s)";
	}
	for (Node::List::Iterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		cobo->setInitDone(done);
	}
}

/**
 * Sets @em mutantMode bit of all CoBo boards to 0 or 1 depending on whether they will be MuTanT driven.
 */
void SystemManager::setCoBoMutantMode(bool enabled)
{
	if (enabled)
	{
		LOG_INFO() << "Setting 'mutantMode' bit of CoBo board(s)";
	}
	for (Node::List::Iterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		cobo->setMutantMode(enabled);
	}
}

/**
 * Sets @em l2Trigger bit of all CoBo boards to 0 or 1 depending in whether they will be MuTanT driven.
 */
void SystemManager::setCoBoLevel2Trigger(CCfg::CConfig globalCfg)
{
	// Enable Level 2 trigger in CoBo boards if enabled in master MuTanT module.
	bool enabled = false;
	if (hasMutantMaster())
	{
		enabled = mutantMaster().nodeCfg(globalCfg)("Level", 2)("isEnabled").getBoolValue(false);
	}

	if (enabled)
	{
		LOG_INFO() << "Setting 'l2trigger' bit of CoBo board(s)";
	}

	for (Node::List::Iterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		cobo->setLevel2Trigger(enabled);
	}
}

/**
 * Looks up which CoBo boards user wants to use.
 */
void SystemManager::selectCoBoBoards(CCfg::CConfig globalCfg)
{
	for (Node::List::Iterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		cobo->setEnabled(cobo->nodeCfg(globalCfg)("isActive").getBoolValue());
	}
}

/**
 * Inter-board communications, clock synchronizations, calibrations.
 * @param config Configuration to use.
 * Prepares all nodes:
 *   - prepare each CoBo board
 *   - prepares each MuTanT module
 *   - find out which slots the CoBo boards are in, using Mutant's alignment procedures
 *   - run MuTanT alignment procedures
 */
void SystemManager::prepare(CCfg::CConfig globalCfg)
{
	// Store number of MuTanT nodes in every MuTanT node
	for (Node::List::Iterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		mutant->setModuleCount(getList<MutantNode>().size());
	}
	// Prepare MuTanT modules
	getList<MutantNode>().doLoop_prepare(globalCfg);

	// Check in configuration which CoBo boards user wants to actually use
	selectCoBoBoards(globalCfg);

	// Identify AMC slots CoBo boards are in
	rc::AmcSlotManager slotManager(*this);
	slotManager.findOccupiedSlots(amcSlots);

	// Run alignment procedures inside each crate
	slotManager.measureDelays();

	// Run alignment procedure between crates
	alignDelaysBetweenChassis();

	// Set all MuTanT modules back to their "normal" mode (CLOCK_TREE_CONFIG2<3>=1)
	setWscaClk2(true);

	// Prepare CoBo boards
	getList<CoBoNode>().doLoop_prepare(globalCfg);
}

/**
 * Sets CLOCK_TREE_CONFIG2<3> in all MuTanT modules.
 *
 * G.Wittwer, Dec 8 2016: This bit must be set to 1 in all cases (master or slave crate) after alignment (mono or multi-crate),
 *  otherwise this fundamental signal (sampling if 1, dead time if 0) is never transfered to CoBo's.
 * G. Wittwer, Dec 13 2016: It should be set to 1 only after both the alignment procedure AND the synchronization phase of CoBo PLL's,
 * otherwise this SYNC_PULSE is not transfered to the CoBo backplane.
 */
void SystemManager::setWscaClk2(bool enabled)
{
	for (Node::List::Iterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		mutant->setWscaClk2(enabled);
	}
}

/**
 * Trigger SYNC_PULSE procedure to be sure that all PLL output frequencies are also phase aligned for all the CoBo's.
 *
 * From G. Wittwer (Dec 8, 2016):
 * Once the alignment procedure is finished, with bit CLOCK_TREE_CONFIG2<3> set back to 1, the CoBo's must not be started,
 * but bit CLOCK_TREE_CONFIG2<9> known as START_PULSE is also used in this case as SYNC_PULSE (see bottom of page 10 in Slow Control Registers Manual).
 * The PLL's of the CoBo's that receive GMC and provide different output frequencies must be aligned in phase for each frequency
 * (e.g. the 25 MHz of CoBo 'n' with the 25 MHz of CoBo 'n+1', etc...).
 * This SYNC_PULSE bit must then be set back to 0, before starting the CoBo's before starting the triggers.
 *
 * This procedure is necessary after alignment and before starting the system. On the CoBo side, during this procedure,
 * register 'mutantConfig' must have its 'moutantMode' bit set to 1 but 'InitDone' should be left to 0 so that CoBo can distinguish this brief passage
 * to 1 of WSCA from a start of the acquisition.
 *
 * It should be done for all MuTanT modules, the order does not matter but it should happen:
 *  - after the alignment phases
 *  - when all GMC and WSCA outputs (register CLOCK_TREE_CONFIG1) are enabled in MuTanT.
 *
 * From G. Wittwer (Dec 13, 2016):
 * For the SYNC_PULSE to be transmitted to CoBo backplanes, CLOCK_TREE_CONFIG2<3> must be set to 0. But it should be set back to 1 afterwards so that CLK2 is distributed to enabled AMC's.
 *
 */
void SystemManager::sendSyncPulses()
{
	if (getList<MutantNode>().size() <= 0)
		return;
	LOG_DEBUG() << "Sending SYNC_PULSE to align phases of all CoBo PLL's...";
	// Set CLOCK_TREE_CONFIG2<9> back to 1
	for (Node::List::Iterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		mutant->restoreClockShelfMasks();
		mutant->setSyncPulse(true);
	}

	// Set CLOCK_TREE_CONFIG2<9> back to 0
	for (Node::List::Iterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		mutant->setSyncPulse(false);
	}
}

/**
 * Run MuTant alignment procedure between multiple chassis.
 *
 * Alignment inside a single chassis is done with a precision of ~300 ps.
 * Alignment between two chassis is done within ~700 ps.
 *
 * After triggering the internal alignment procedure performed on each chassis, the ECC should, if there is more than 1 MuTanT module:
 *   - Set all active CoBos? boards into "Normal" MuTanT-link mode.
 *   - For each MuTanT slave:
 *       - For each clock (GMC, WSCA):
 *          - master.step1(clock)
 *          - slave.step1(clock)
 *          - master.step2(clock)
 *          - slave.step2(clock)
 *          - master.step3(clock)
 *          - slave.step3(clock)
 *       - Print alignment precision?
 *   - Switch all MuTanT to their normal state: CLOCK_TREE_CONFIG2<3>=1
 *   - Switch all CoBo boards to their normal MuTanT-link mode.
 *
 * Before alignment, in the prepare step, master and slaves should already be partially configured:
 *   - Set master/slave bits, VPA_ON, VPB_ON in GLOBAL_CONF
 *   - Set GMC_WSCA_SEL/PromptWSCA to 0 in slaves, to 1 in master if there is at least 1 slave, to 0 in master otherwise
 */
void SystemManager::alignDelaysBetweenChassis()
{
	if (getList<MutantNode>().size() <= 1u)
		return;

	LOG_DEBUG() << "Aligning delays between shelves...";	
	rc::AmcSlotManager slotManager(*this);
	try
	{
		// Set all active CoBo boards in loopback mode, others in Disconnected mode
		// That way, these 'inactive' boards will not appear in CLOCK_TREE_CONFIG1/2 registers
		slotManager.setCoBoMutantLink(cobo::Disconnected, false);
		slotManager.setCoBoMutantLink(cobo::Normal, true);

		// Loop over slaves
		MutantSeq slaves(mutantSlaves());
		for (size_t i=0; i < slaves.size(); ++i)
		{
			// Loop over two clocks
			for (size_t clock_num = 0; clock_num < 2; ++clock_num)
			{
				// Set CoBo-MuTanT link mode
				cobo::MutantLinkMode mode = static_cast< cobo::MutantLinkMode >(cobo::GMC + clock_num);
				slotManager.setCoBoMutantLink(mode, true);

				// Loop over three steps
				mutant::Clock c = static_cast< mutant::Clock >(mutant::GMC + clock_num);
				for (size_t step = 0; step < 3; ++step)
				{
					alignDelaysBetweenChassis(*slaves[i], c, step);
				}
			}
		}

		// Set all MuTanT modules back to their "normal" mode (CLOCK_TREE_CONFIG2<3>=1)
		setWscaClk2(true);

		// Set all CoBo boards in normal mode
		slotManager.setCoBoMutantLink(cobo::Normal, false);

	}
	catch (const mutant::AlignmentAlarm & e)
	{
		// Set all MuTanT modules back to their "normal" mode (CLOCK_TREE_CONFIG2<3>=1)
		setWscaClk2(true);
		// Set all CoBo boards in normal mode
		slotManager.setCoBoMutantLink(cobo::Normal, false);
		throw CmdException(e.reason);
	}
	catch (const CmdException & e)
	{
		// Set all MuTanT modules back to their "normal" mode (CLOCK_TREE_CONFIG2<3>=1)
		setWscaClk2(true);
		// Set all CoBo boards in normal mode
		slotManager.setCoBoMutantLink(cobo::Normal, false);
		throw;
	}
}
//__________________________________________________________________________________________________
void SystemManager::alignDelaysBetweenChassis(MutantNode & slave, get::mutant::Clock c, size_t step)
{
	LOG_DEBUG() << "Alignment step no. " << step << " for MuTanT module " << slave.moduleId_;
	mutant::MultiShelveAlignmentResults r;
	bool verbose = true;

	if (step <= 0)
	{
		// Run step on master
		mutantMaster().alignmentServer()->multiShelveAlignDelayStep1(c, verbose, r);
		if (not r.comment.empty())
		{
			LOG_INFO() << r.comment;
		}

		// Run step on slave
		slave.alignmentServer()->multiShelveAlignDelayStep1(c, verbose, r);
		if (not r.comment.empty())
		{
			LOG_INFO() << r.comment;
		}
	}
	else if (step == 1)
	{
		// Run step on master
		mutantMaster().alignmentServer()->multiShelveAlignDelayStep2(c, verbose, r);
		if (not r.comment.empty())
		{
			LOG_INFO() << r.comment;
		}

		// Run step on slave
		slave.alignmentServer()->multiShelveAlignDelayStep2(c, verbose, r);
		if (not r.comment.empty())
		{
			LOG_INFO() << r.comment;
		}
	}
	else
	{
		// Run step on master
		mutantMaster().alignmentServer()->multiShelveAlignDelayStep3(c, verbose, r);
		if (not r.comment.empty())
		{
			LOG_INFO() << r.comment;
		}

		// Run step on slave
		slave.alignmentServer()->multiShelveAlignDelayStep3(c, verbose, r);
		if (not r.comment.empty())
		{
			LOG_INFO() << r.comment;
		}
	}
}

/**
 * Powers off all AsAd boards of all CoBo boards.
 */
void SystemManager::unprepare()
{
	CCfg::Object dummyAttribute("DummyConfiguration");
	getList<CoBoNode>().doLoop_unprepare(CCfg::CConfig(&dummyAttribute)); // With empty config
}

/**
 * Loads physics parameters into the registers of the boards (thresholds, ...).
 */
void SystemManager::configure(CCfg::CConfig& globalCfg, const std::string & dataLinksString)
{
	getList<CoBoNode>().doLoop_configure(globalCfg);

	// Enable Level 2 trigger in CoBo boards if enabled in master MuTanT module.
	setCoBoLevel2Trigger(globalCfg);

	// Tell CoBo boards whether they will be MuTanT driven
	// BR. 2014/10/08  action displaced at the end of configure, because problem with ASAD config if this bit is set
	setCoBoMutantMode(getList<MutantNode>().size() > 0);

	// Set CLOCK_TREE_CONFIG2<3> to 0 so that SYNC_PULSE will be sent instead of CLK2
	setWscaClk2(false);
	// Send SYNC_PULSE signals to align CoBo PLL's
	sendSyncPulses();
	// Set CLOCK_TREE_CONFIG2<3> back to 1 before the start of acquisition
	setWscaClk2(true);

	// Make all active CoBo boards send 'ready' status to MuTanT modules.
	// Must be done before configuring MuTanT modules as it affects MuTanT B's MRDY status bit.
	setCoBoInitDone(true);

	// Configure all MuTanT modules
	getList<MutantNode>().doLoop_configure(globalCfg);

	// Check coherence between configurations of the 3 MuTanT modules
	// TODO: check s0rdy and s1rdy here ?
	checkCoherenceBetweenMutantModules();

	// Check coherence between configuration of CoBo boards and MuTanT modules
	checkCoBoMutantCoherence();

	// TODO: disable AMC slots if CoBo not active?

	// Parse list of data links
	DataLinkSet dataLinks;
	try
	{
		dataLinks.fromXml(dataLinksString);
	}
	catch (const std::exception & e)
	{
		throw CmdException(std::string("Error parsing data links '") + dataLinksString + "': " + e.what());
	}

	// Set data links for CoBo boards
	LOG_INFO() << "Creating links between boards and data routers...";
	LOG_DEBUG() << "Available data links:\n" << dataLinks;

	for (Node::List::Iterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		cobo->setDataLinks(dataLinks); // Set data router parameters
		cobo->daqConnect(); // Connect board to data router
	}

	// BR 2015/03 Set data links for MuTanT modules
	for (Node::List::Iterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		if (mutant->isReadoutEnabled())
		{
			mutant->setDataLinks(dataLinks); // Set data router parameters
			mutant->daqConnect(); // Connect module to data router
		}
	}

	// TODO: DataRouter::runStart() ?
}

/**
 * Breaks any open connection between boards and data routers.
 */
void SystemManager::breakup()
{
	for (Node::List::Iterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		cobo->daqDisconnect(); // Disconnect board from data router
	}

	for (Node::List::Iterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		mutant->daqDisconnect();
	}

	// Make all CoBo boards send <em>not configured</em> status to MuTanT modules.
	setCoBoInitDone(false);

	// Go back to Mutant-less CoBo mode to allow configuration of CoBo
	// Even though "the timing of setting the bit should not matter"...
	setCoBoMutantMode(false);
}

/**
 *  Checks coherence of configurations of the MuTanT modules
 */
void SystemManager::checkCoherenceBetweenMutantModules()
{
	if (getList<MutantNode>().size() > 0)
	{
		// Check number of MuTanT master modules
		size_t masterCount = 0;
		for (Node::List::Iterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
		{
			if (not mutant->isSlaveModule()) masterCount++;
		}

		// At least one master, at most 1 master
		if (masterCount != 1)
		{
			std::ostringstream oss;
			oss << masterCount << " of the " << getList<MutantNode>().size() << " MuTanT module(s) have been declared as master!";
			throw CmdException(oss.str());
		}

		// Check number of slaves
		MutantSeq slaves(mutantSlaves());
		// At most 2 slaves
		if (slaves.size() > 2)
		{
			throw CmdException("There can be at most 2 MuTanT slave modules!");
		}
		else if (slaves.size() == 2)
		{
			// Slaves must have different slave identifiers
			if (slaves[0]->config.slaveId == slaves[1]->config.slaveId)
			{
				throw CmdException("MuTanT slave modules cannot share the same slave identifier!");
			}
			// Slaves must be connected to different ports
			if (slaves[0]->config.isSlavePortA and slaves[1]->config.isSlavePortA)
			{
				throw CmdException(std::string("MuTanT slave modules cannot be both connected to port ") + (slaves[0]->config.isSlavePortA?"A!":"B!"));
			}
		}

		// Check that slave has SNUM=0 if in port A, SNUM=1 if in port B
		for (size_t i=0; i < slaves.size(); ++i)
		{
			const MutantConfigParser & cfg = slaves[0]->config;
			if ((cfg.isSlavePortA and cfg.slaveId == 1) or (cfg.isSlavePortA == false and cfg.slaveId == 0))
			{
				std::ostringstream oss;
				oss << "Slave with number " << cfg.slaveId << " cannot be assigned to port " << (cfg.isSlavePortA?"A":"B") << "!";
				throw CmdException(oss.str());
			}
		}
	}
}

/**
 * Checks CoBo boards and MuTanT modules are configured coherently.
 */
void SystemManager::checkCoBoMutantCoherence()
{
	// TODO: coherence between CoBo and MuTanT config (sourceId...)
}

void SystemManager::start()
{
	LOG_INFO() << "Starting new run";
	// Loop over CoBo boards and run pre-start actions
	for (Node::List::Iterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		cobo->runInit();
	}

	// Start CoBo boards
	CCfg::Object dummyAttribute("DummyConfiguration");
	getList<CoBoNode>().doLoop_start(CCfg::CConfig(&dummyAttribute)); // empty config for start.

	if (getList<MutantNode>().size() > 0)
	{
		// Loop over slave MuTanT modules
		MutantSeq slaves(mutantSlaves());
		for (size_t i=0; i < slaves.size(); ++i)
		{
			LOG_DEBUG() << "Starting MuTanT slave module " << slaves[i]->moduleId();
			// Clear counters and scalers before start
			slaves[i]->runInit();
			// Start slave
			slaves[i]->impl_start();
		}
		// Clear master's counters and scalers before start
		LOG_DEBUG() << "Starting MuTanT master module " << mutantMaster().moduleId();
		mutantMaster().runInit();
		// Start master MuTanT module
		mutantMaster().impl_start();
	}
}

void SystemManager::stop()
{
	LOG_INFO() << "Stopping run";

	// B.R. 2014/10/10 Stop Cobo first until upgrade of the MUTANT firmware
	// G.W. 2016/03/16 Stop MuTanT first to avoid "shelf is not ready!" errors
	if (getList<MutantNode>().size() > 0)
	{
		// Stop master MuTanT module
		mutantMaster().impl_stop();

		// Stop slave MuTanT modules
		MutantSeq slaves(mutantSlaves());
		for (size_t i=0; i < slaves.size(); ++i)
		{
			slaves[i]->impl_stop();
		}
	}

	// Stop CoBo boards
	CCfg::Object dummyAttribute("DummyConfiguration");
	getList<CoBoNode>().doLoop_stop(CCfg::CConfig(&dummyAttribute)); // With empty CCfg
}

void SystemManager::pause()
{
	LOG_INFO() << "Pausing run";
	if (getList<MutantNode>().size() > 0)
	{
		// Stop master MuTanT module
		mutantMaster().impl_stop();

		// Stop slave MuTanT modules
		MutantSeq slaves(mutantSlaves());
		for (size_t i=0; i < slaves.size(); ++i)
		{
			slaves[i]->impl_stop();
		}
	}

	// Stop CoBo boards
	CCfg::Object dummyAttribute("DummyConfiguration");
	getList<CoBoNode>().doLoop_stop(CCfg::CConfig(&dummyAttribute)); // With empty CCfg
}

void SystemManager::resume()
{
	LOG_INFO() << "Resuming run";
	// Resume CoBo boards
	getList<CoBoNode>().doLoop_start(CCfg::CConfig());

	if (getList<MutantNode>().size() > 0)
	{
		// Resume slave MuTanT modules
		MutantSeq slaves(mutantSlaves());
		for (size_t i=0; i < slaves.size(); ++i)
		{
			slaves[i]->impl_start();
		}
		// Resume master MuTanT module
		mutantMaster().impl_start();
	}
}

/**
 * Returns the list of data senders.
 */
DataSenderSet SystemManager::dataSenders() const
{
	DataSenderSet senders;
	// Loop over CoBo boards

	for (Node::List::ConstIterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		// Add board links
		const DataSenderSet boardSenders = cobo->dataSenders();
		senders.insert(boardSenders.begin(), boardSenders.end());
	}
	// Loop over MuTanT modules
	for (Node::List::ConstIterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		// Add board links
		const DataSenderSet boardSenders = mutant->dataSenders();
		senders.insert(boardSenders.begin(), boardSenders.end());
	}
	return senders;
}

/**
 * Returns the mapping between CoBo or MuTanT boards and their data routers.
 */
DataLinkSet SystemManager::dataLinks() const
{
	DataLinkSet links;
	// Loop over CoBo boards
	for (Node::List::ConstIterator<CoBoNode> cobo = iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		// Add board links
		const DataLinkSet boardLinks = cobo->dataLinks();
		links.insert(boardLinks.begin(), boardLinks.end());
	}
	// Loop over MuTanT modules
	for (Node::List::ConstIterator<MutantNode> mutant = iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		// Add board links
		const DataLinkSet boardLinks = mutant->dataLinks();
		links.insert(boardLinks.begin(), boardLinks.end());
	}
	return links;
}


} /* namespace rc */
} /* namespace get */
