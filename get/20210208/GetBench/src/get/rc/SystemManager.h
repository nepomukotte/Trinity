/**
 * @file SystemManager.h
 * @date Mar 18, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: SystemManager.h 1738 2017-11-23 15:50:07Z psizun $
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

#ifndef get_rc_SystemManager_h_INCLUDED
#define get_rc_SystemManager_h_INCLUDED

#include "get/rc/AmcSlotManager.h"
#include "get/rc/ConfigId.h"
#include "get/rc/DataSenderSet.h"
#include "mdaq/rc/HardwareNode.h"
#include <CCfg/CConfig.h>
#include <boost/ptr_container/ptr_map.hpp>
#include <rc/CoBoNode.h>
#include <rc/MutantNode.h>
#include <map>

namespace get
{
namespace rc
{
class BackEnd;
//_________________________________________________________________________________________________
/**
 * A class to manage a full GET system: a set of CoBo board(s) and Mutant module(s).
 */
class SystemManager : public Node::ListList
{
public:
	SystemManager(BackEnd& backend);
	~SystemManager();

	void removeNodes();
	CoBoNode & getCoBoNode(const std::string & name);

	/// @name State machine transitions
	/// @{
	void describe();
	void prepare(CCfg::CConfig globalCfg);
	void unprepare();
	void configure(CCfg::CConfig& globalCfg, const std::string& dataLinksString);
	void breakup();
	void start();
	void stop();
	void pause();
	void resume();
	/// @}

	void setWscaClk2(bool enabled);
	void sendSyncPulses();
	void alignDelaysBetweenChassis();
	void alignDelaysBetweenChassis(MutantNode & slave, get::mutant::Clock c, size_t step);
	MutantNode & mutantMaster();
	bool hasMutantMaster();
	MutantSeq mutantSlaves();
	void checkCoherenceBetweenMutantModules();

	void selectCoBoBoards(CCfg::CConfig config);
	void setCoBoInitDone(bool done);
	void setCoBoMutantMode(bool enabled);
	void setCoBoMutantLink(const cobo::MutantLinkMode & mode, bool skipDisabledBoards);
	void setCoBoLevel2Trigger(CCfg::CConfig config);
	void checkCoBoMutantCoherence();

	DataSenderSet dataSenders() const;
	DataLinkSet dataLinks() const;

	friend class ::get::rc::AmcSlotManager;
protected:
	std::map< std::string, AmcSlotManager::SlotUID > amcSlots; ///< Map of AMC slots occupied by CoBo boards.
};
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
#endif /* get_rc_SystemManager_h_INCLUDED */
