/**
 * @file MutantNode.h
 * @date Sep 18, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: MutantNode.h 1799 2021-02-08 15:42:53Z psizun $
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

#ifndef get_rc_MutantNode_h_INCLUDED
#define get_rc_MutantNode_h_INCLUDED

#include "MutantConfigParser.h"
#include "get/mutant/AlignmentServer.h"
#include <boost/ptr_container/ptr_map.hpp>
#include <rc/Node.h>
#include <map>
#include <string>
#include <vector>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
/**
 * A class to manage a Mutant module.
 *
 * @see Mutant Slow Control Registers (v1.4)
 *
 */
class MutantNode : public Node
{
public:
	static const std::string NODE_TYPE;
	typedef boost::ptr_map<std::string, MutantNode> Map;

public:
	MutantNode(const std::string & moduleId="", mdaq::rc::HardwareNode* node=0, get::rc::BackEnd* backend=0);
	std::string moduleId() const { return id().instanceId; }
	void runInit();
	virtual DataSenderSet dataSenders() const;
	bool isSlaveModule() const;
	bool isReadoutEnabled() const;
	bool isTriggerEnabled() const;
	void setModuleCount(size_t n);

protected:
	virtual void impl_prepare(CCfg::CConfig nodeCfg);
	virtual void impl_unprepare();
	virtual void impl_configure(CCfg::CConfig nodeCfg);
	virtual void impl_start();
	virtual void impl_stop();

protected:
	/// @name Methods for preparing MuTanT module
	/// @{
	void checkPowerSupply();
	void checkIdentifier();
	void ensureModuleNotRunning();
	/// @}
	/// @name Methods for configuring MuTanT module
	/// @{
	void configureRegistersMutantA(CCfg::CConfig cfg);
	void setWscaClk2(bool enabled);
	void setSyncPulse(uint8_t sync_pulse);
	void configureTdcStatAndCardBFlags(CCfg::CConfig cfg);
	void configureRegistersMutantB(CCfg::CConfig cfg);
	void registerBeastBoards(CCfg::CConfig /* cfg */);
	void tuneVcxoFrequency(const uint8_t vcxo_control = 0x80);
	void selectClockSource(const bool & isSlave, const bool & slavePortA);
	void setGlobalConfigurationRegister();
	void selectReadoutMode();
	void selectTestMode();
	void setTimestampAndEventCounters();
	void clearEventCounters();
	void enableScaler(const size_t & scalerId, const bool enable);
	void enableScalers(bool enable);
	void clearScalers();
	void configureScalers();
	void setSourceId(const uint8_t id);
	void setFrameType(const uint8_t type);
	void storeClockShelfMasks();
	void setClockShelfMasks(const uint16_t & gmcMask, const uint16_t & wscaMask);
	void restoreClockShelfMasks();
	void setShelfMask();
	void setAmcTagOnly();
	void setUserRegister();
	void configureLevel0Trigger();
	void configureLevel1Trigger();
	void configureLogicalInspectionSignals();
	void checkStatus();
	void checkAmcStatus();
	void setModuleRunning();
	/// @}
	const MutantConfigParser & configData() { return config; };
	DataSenderId dataSenderId() const;
	mutant::AlignmentServerPrx alignmentServer();
	friend class AmcSlotManager;
	friend class SystemManager;
protected:
	std::string moduleId_; ///< Mutant module identifier.
	MutantConfigParser config; ///< Parsed configuration.
	size_t mutantModuleCount_; ///< Number of MuTanT modules.
	uint16_t fMaskGmc_;
	uint16_t fMaskWsca_;
};
typedef std::vector< MutantNode* > MutantSeq;
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
#endif /* get_rc_MutantNode_h_INCLUDED */
