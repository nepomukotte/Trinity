/**
 * @file AmcSlotManager.cpp
 * @date Feb 24, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: AmcSlotManager.cpp 1727 2017-11-14 12:34:14Z psizun $
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

#include "AmcSlotManager.h"
#include "get/cobo/MutantLinkManager.h"
#include "utl/Logging.h"
#include <cmath>
#include <iomanip>
#include <utility>
#include <vector>
using mdaq::utl::CmdException;

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
AmcSlotManager::AmcSlotManager(Node::ListList& nodeLists) : nodeLists(nodeLists)
{
}
//_________________________________________________________________________________________________
AmcSlotManager::SlotUID AmcSlotManager::buildSlotUID(const NodeId & id, const size_t & slot)
{
	return std::make_pair(id, slot);
}
//_________________________________________________________________________________________________
std::string AmcSlotManager::buildSlotUIDString(const NodeId & id, const size_t & slot)
{
	std::ostringstream oss;
	oss << id.toString() << ".Slot[" << std::setw(2) <<  (1 + slot) << ']';
	return oss.str();
}
//_________________________________________________________________________________________________
/**
 * Sets all CoBo boards in the same Mutant link mode.
 * @param mode Mutant link mode to set the CoBo boards in.
 * @param skipDisabledBoards Whether to skip boards previously disabled.
 */
void AmcSlotManager::setCoBoMutantLink(const cobo::MutantLinkMode & mode, bool skipDisabledBoards)
{
	for (Node::List::Iterator<CoBoNode> cobo = nodeLists.iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		// Skip CoBo boards where isActive was set to false in configuration file
		if (skipDisabledBoards and not cobo->isEnabled())
			continue;
		cobo->mutantLinkManager()->setupMutantLink(mode);
	}
}
//_________________________________________________________________________________________________
/**
 * Finds AMC slot occupied by CoBo board with given ID by triggering the alignment procedure of all
 * MuTanT modules after setting the coBo board in given feedback mode.
 *
 * All boards should be in normal Mutant link mode prior and after this search.
 *
 * @param coboId Identifier of the CoBo board to localize
 * @param clock Feedback mode to use for the search.
 */
AmcSlotManager::SlotUID AmcSlotManager::findSlotOfCoBo(const std::string & coboId, const mutant::Clock & clock)
{
	// Enter feedback mode
	cobo::MutantLinkMode mode = static_cast< cobo::MutantLinkMode >(cobo::GMC + clock);
	nodeLists.getNode<CoBoNode>(coboId).mutantLinkManager()->setupMutantLink(mode);

	// Loop over MuTanT modules
	std::map< SlotUID, mutant::SlotAlignmentResult > results;
	for (Node::List::Iterator<MutantNode> mutant = nodeLists.iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		// Run alignment procedure
		mutant::AlignmentResults r;
		mutant->alignmentServer()->alignDelayTDC(clock, false, r);

		// Store results where slot has been found busy
		mutant::AlignmentResults::const_iterator resultIter;
		for (resultIter = r.begin(); resultIter != r.end(); ++resultIter)
		{
			if ((*resultIter).second.isAmcPresent)
			{
				results[buildSlotUID(mutant->id(), 1 + resultIter->first)] = (*resultIter).second;
			}
		}
	}

	// Restore normal mode
	nodeLists.getNode<CoBoNode>(coboId).mutantLinkManager()->setupMutantLink(cobo::Normal);

	// Process results
	if (results.size() != 1)
	{
		std::ostringstream msg;
		msg << "Could not identify AMC slot occupied by CoBo[" << coboId << "]. Found " << results.size() << " slot(s) in clock mode " << clock;
		throw CmdException(msg.str());
	}

	return results.begin()->first;
}
//_________________________________________________________________________________________________
/**
 * Identifies slots occupied by CoBo boards, both in GMC and WSCA feedback modes.
 */
void AmcSlotManager::findOccupiedSlots(std::map< std::string, SlotUID > & slots)
{
	if (nodeLists.getList<MutantNode>().size() <= 0)
	{
		LOG_DEBUG() << "Skipping localization of CoBo board(s) since no MuTanT module is available";
		return;
	}
	LOG_INFO() << "Searching for the slots occupied by CoBo boards";

	// Set all boards in normal mode
	setCoBoMutantLink(cobo::Normal, false);

	slots.clear();
	for (Node::List::Iterator<CoBoNode> cobo = nodeLists.iterate<CoBoNode>(); cobo.hasNext(); ++cobo)
	{
		// Skip CoBo boards where isActive was set to false in configuration file
		if (not cobo->isEnabled()) continue;

		//LOG_DEBUG() << "Identifying AMC slot of CoBo board '" << cobo.id().toString() << "'";

		// Find slot using GMC feedback mode
		SlotUID gmcSlot = findSlotOfCoBo(cobo->id().instanceId, mutant::GMC);

		// Find slot using WSCA feedback mode
		SlotUID wscaSlot = findSlotOfCoBo(cobo->id().instanceId, mutant::WSCA);

		// Check we found the same slot with each feedback mode
		if (gmcSlot != wscaSlot)
		{
			std::ostringstream msg;
			msg << cobo->id().toString() << " was localized "
				<< " in slot " << gmcSlot.second  << " of " << gmcSlot.first.toString()  << "'s crate using GMC feedback mode and "
				<< " in slot " << wscaSlot.second << " of " << wscaSlot.first.toString() << "'s crate using  MC feedback mode";
			throw CmdException(msg.str());
		}

		LOG_DEBUG() << cobo->id().toString() << " is located in AMC slot " << gmcSlot.second
			        << " of " << gmcSlot.first.toString() << "'s crate";
		slots[cobo->id().instanceId] = gmcSlot;
	}
}
//_________________________________________________________________________________________________
/**
 * Runs alignment procedures in both loopback modes with all CoBo boards and MuTanT modules.
 */
void AmcSlotManager::measureDelays()
{
	if (nodeLists.getList<MutantNode>().size() <= 0)
	{
		LOG_DEBUG() << "Skipping alignment procedures since no MuTanT module is available";
		return;
	}

	// Run alignment procedure in GMC loopback mode
	AlignmentResults gmcResults;
	measureDelays(mutant::GMC, gmcResults);

	// Run alignment procedure in WSCA loopback mode
	AlignmentResults wscaResults;
	measureDelays(mutant::WSCA, wscaResults);

	// Store masks of AMC slots found occupied
	for (Node::List::Iterator<MutantNode> mutant = nodeLists.iterate<MutantNode>(); mutant.hasNext(); ++mutant)
	{
		mutant->storeClockShelfMasks();
	}

	// Check coherence of results
	if (gmcResults.size() != wscaResults.size())
		throw CmdException("GMC and WSCA alignment procedures did not return the same number of results!");
	AlignmentResults::const_iterator gmcIt;
	for (gmcIt = gmcResults.begin(); gmcIt != gmcResults.end(); ++gmcIt)
	{
		if (wscaResults.find(gmcIt->first) == wscaResults.end())
		{
			std::ostringstream oss;
			oss << "AMC slot " << buildSlotUIDString(gmcIt->first.first, gmcIt->first.second) << " was found in GMC mode but not in WSCA mode!";
			throw CmdException(oss.str());
		}
		else if (wscaResults[gmcIt->first].isAmcPresent != gmcIt->second.isAmcPresent)
		{
			std::ostringstream oss;
			oss << "AMC slot " << buildSlotUIDString(gmcIt->first.first, gmcIt->first.second) << " is "
					<< (gmcIt->second.isAmcPresent?"present":"absent") << " according to GMC mode but not to WSCA mode!";
			throw CmdException(oss.str());
		}
		//Check delays do not differ by more than 10 ns; the amcDeltaDelay delay value is expressed in pico-seconds
		else if (gmcIt->second.isAmcPresent and std::fabs(gmcIt->second.amcDeltaDelay - wscaResults[gmcIt->first].amcDeltaDelay) > 10000)
		{
			throw CmdException(std::string("GMC and WSCA delay values for slot ") + buildSlotUIDString(gmcIt->first.first, gmcIt->first.second) + " differ by more than 10ns!");
		}
	}

}
//_________________________________________________________________________________________________
/**
 * Runs alignment procedure in given loopback mode with all Cobo boards and MuTanT modules.
 * @param clock Loopback mode.
 * @param results Results of alignment procedures
 */
void AmcSlotManager::measureDelays(const mutant::Clock & clock, AlignmentResults & results)
{
	LOG_INFO() << "Running Mutant alignment procedure in mode " << (clock == mutant::GMC?"GMC":"WSCA");

	// Set all active CoBo boards in loopback mode, others in Disconnected mode
	// That way, these 'inactive' boards will not appear in CLOCK_TREE_CONFIG1/2 registers
	setCoBoMutantLink(cobo::Disconnected, false);
	cobo::MutantLinkMode mode = static_cast< cobo::MutantLinkMode >(cobo::GMC + clock);
	setCoBoMutantLink(mode, true);

	// Loop over MuTanT modules
	try
	{
		results.clear();
		for (Node::List::Iterator<MutantNode> mutant = nodeLists.iterate<MutantNode>(); mutant.hasNext(); ++mutant)
		{
			// Run alignment procedure
			mutant::AlignmentResults r;
			mutant->alignmentServer()->alignDelayTDC(clock, false, r);

			// Copy and log results
			mutant::AlignmentResults::const_iterator resultIter;
			for (resultIter = r.begin(); resultIter != r.end(); ++resultIter)
			{
				results[buildSlotUID(mutant->id(), resultIter->first)] = resultIter->second;

				LOG_DEBUG() << "Clock[" << (clock == mutant::GMC?"GMC":"WSCA") << "] "
					<< buildSlotUIDString(mutant->id(), resultIter->first)
					<< " isAmcPresent=" << (*resultIter).second.isAmcPresent
					<< " deltaDelay=" << (*resultIter).second.amcDeltaDelay
					<< " progDelay=" << (*resultIter).second.amcProgDelay;
			}
		}
	}
	catch (const CmdException & e)
	{
		// Set all CoBo boards in normal mode
		setCoBoMutantLink(cobo::Normal, false);
		throw;
	}

	// Set all CoBo boards in normal mode
	setCoBoMutantLink(cobo::Normal, false);
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */

