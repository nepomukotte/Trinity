/* =====================================================================================================================
 * BackEnd.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class get::rc::BackEnd
 * Created on: Apr 4, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GET data acquisition framework.
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
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef get_rc_BackEnd_h_INCLUDED
#define get_rc_BackEnd_h_INCLUDED

#include "get/DefaultPortNums.h"
#include "get/GetEcc.h"
#include "get/rc/SM.h"
#include "get/rc/ConfigId.h"
#include "get/rc/SystemManager.h"
#include "get/mt/AlarmService.h"
#include <mdaq/Ecc.h>
#include <mdaq/EccBackEnd.h>
#include <dhsm/Engine.h>
#include <utl/BitFieldRegister.hpp>
#include <utl/net/SocketAddress.h>
#include ICE_HEADER(Communicator.h)
#include <boost/cstdint.hpp>
#include <set>
#include <string>

namespace get
{
namespace rc
{
//__________________________________________________________________________________________________
// Forward declarations
	class ConfigManager;
	class ConfigIdSet;
	class SubConfigIdSet;
//__________________________________________________________________________________________________
/**
 * Backend of the GET ECC server.
 */
class BackEnd : public virtual ::mdaq::EccBackEnd, protected dhsm::StateData
{
public:
	BackEnd();
	virtual ~BackEnd();
	void setConfigManager(std::auto_ptr< ConfigManager > manager);
	void setConfigManager(const std::string & repoUrl = ".");
	const ConfigManager & configManager() const;
	void daqStart();
	void daqStop();

	/// @name Alarms
	/// @{
	cobo::AsAdPulserMgrPrx asadPulserMgr();
	cobo::AsAdAlarmMonitorPrx asadAlarmMonitor();
	void createAlarmLogger(const ::utl::net::SocketAddress & address = ::utl::net::SocketAddress(Default::alarmLoggerPortNum));
	mt::AlarmCallbackPrx alarmLogger();
	mt::AlarmServicePrx alarmService();
	void subscribeAlarms(bool enabled);
	/// @}

	cobo::LedManagerPrx ledManager();
	void asadSetISPM(const std::string & moduleName, const std::string & isp1SignalName, const std::string & isp2SignalName);
	void connect(const std::string& objectType, const std::string& objectName, bool isBatched=false);
	virtual void loadHwDescription(const std::string& cfgUrl);
	void getStateMachineStatus(SM::Status& response);
	DataSenderSet dataSenders() const;
	DataLinkSet dataLinks() const;
	std::string specialCmd(const std::string & command);
	/// @name Mutant related methods
	/// @{
	void setupCoBoMutantLink(get::cobo::MutantLinkMode mode);
	mutant::AlignmentServerPrx mutantAlignmentServer();
	/// @}
	/// @name CoBo specific methods
	/// @{
	void coboSetWritingClockFrequency(const float & frequency_MHz, bool pllConfigEnabled=true, const cobo::PllRegisterMap & regs = cobo::PllRegisterMap());
	void coboSetReadingClockFrequency(const float &  frequency_MHz);
	void coboSwitchToSlowControlMode(bool slowMode=true);
	void coboSelectSingleReadoutChannel(const size_t & chanId);
	void coboSelectSingleTestChannel(const size_t & chanId);
	void coboSetAgetChannelLsbThreshold(const size_t & chanId, const uint8_t & threshold);
	/// @}
	/// @name State machine transitions
	/// @{
	void describe(const std::string & subConfigId);
	void prepare(const std::string & subConfigId);
	void reset();
	void configure(const std::string & subConfigId, const std::string &);
	void start();
	void stop();
	void breakup();
	void pause();
	void resume();
	/// @}
	/// @name Management of nodes
	/// @{
	void removeAllHwNodes();
	void describeNode(const std::string nodeType, const std::string & instanceId, const std::string & cfgUrl);
	void prepareNode(const std::string & nodeType, const std::string & instanceId, const std::string & configId);
	void prepareNode(const std::string & nodeType, const std::string & instanceId, const CCfg::CConfig & config);
	void configureNode(const std::string & nodeType, const std::string & instanceId, const std::string & configId);
	void configureNode(const std::string & nodeType, const std::string & instanceId, const CCfg::CConfig & configId);
	/// @}
	/// @name Management of configurations
	/// @{
	void getConfigIds(ConfigIdSet &);
	void getSubConfigIds(const std::string & subConfigType, SubConfigIdSet & subConfigIds) const;
	virtual std::string getConfigSet(const ConfigId & configId) const;
	virtual std::string getLastConfigSet() const;
	virtual std::string addSubConfig(const std::string & subConfigType, const std::string & subConfig, const std::string & subConfigId = "") const;
	virtual std::string getSubConfig(const std::string & subConfigType, const std::string & subConfigId) const;
	/// @}
	/// @name CoBo PLL device
	/// @{
	void loadPLLConfig(const float & ckwFreq_MHz);
	void writePLLRegister(const uint32_t & value);
	void resetPLL();
	/// @}
protected:
	/// @name State machine actions
	/// @{
	void onDescribe(std::string configId);
	void onUnDescribe();
	void onPrepare(std::string configId);
	void onUnPrepare();
	typedef std::pair < std::string, std::string > ConfigParams; ///< configID and dataRouters
	void onConfigure(ConfigParams);
	void onStart();
	void onStop();
	void onBreakup();
	void onPause();
	void onResume();
	/// @}
	std::string buildNodeId(const std::string & nodeType, const std::string & instanceId);
	bool isNodeInstanceActive(CCfg::Object & instance);
private:
	/// @name State machine methods
	/// @{
	SM::State smState();
	inline static SM::State toSmState(const std::string & label);
	SM::Status smStatus();
	void handleStateMachineException(SM::Error errorCode, const std::string & msg);
	void resetErrorFlag();
	void createEngine();
	/// @}
	void resetCoBoReadPointers();
private:
	/// @name State machine attributes
	/// @{
	static const char ST_IDLE[];     ///< Idle state label
	static const char ST_DESCRIBED[];///< Described state label
	static const char ST_PREPARED[]; ///< Prepared state label

	static const char ST_ACTIVE[];   ///< Active state label
	static const char ST_READY[];    ///< Ready state label
	static const char ST_RUNNING[];  ///< Running state label
	static const char ST_PAUSED[];   ///< Paused state label

	static const char EV_PREPARE[];
	static const char EV_DESCRIBE[];
	static const char EV_UNDO[];
	static const char EV_CONFIGURE[];
	static const char EV_START[];
	static const char EV_STOP[];
	static const char EV_PAUSE[];
	static const char EV_RESUME[];
	static const char EV_BREAK[];

	dhsm::Engine engine;
	SM::Error smErrorFlag;
	/// @}
	std::auto_ptr< ConfigManager > configManagerPtr; ///< Accessor to configuration database.
	SystemManager systemMgr_; ///< Experiment setup manager.
	ConfigId lastConfigId; ///< Last configuration identifier.
	::utl::net::SocketAddress alarmLoggerAddress_;
};

} // namespace rc
} // namespace get

#endif /* get_rc_BackEnd_h_INCLUDED */
