/* =====================================================================================================================
 * GetEcc.ice
 * ---------------------------------------------------------------------------------------------------------------------
 * get::GetEcc remote interface definition.
 * Created on: 4 May 2010 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of the Mordicus Real-Time Software
 * Development Framework.
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
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef get_GetEcc_ice_INCLUDED
#define get_GetEcc_ice_INCLUDED

#include "mdaq/Ecc.ice"
#include "get/mt/AlarmService.ice"
#include "get/cobo/MutantLinkManager.ice"
#include "get/cobo/LedManager.ice"
#include "get/mutant/AlignmentServer.ice"

module get
{
		module cobo
		{
			dictionary< short, int > PllRegisterMap;
		};

		/** Run Control **/
		module rc
		{
			enum State { Off, Idle, Described, Prepared, Ready, Running, Paused };
	
			enum Error
			{
				NoErr, WhenDescribe, WhenPrepare, WhenConfigure, WhenStart,
				WhenStop, WhenPause, WhenResume, WhenBreakup, WhenReset
			};
	
			struct Status
			{
				State s;
				Error e;
			};
			
			/** ECC State Machine **/
			interface StateMachine
			{
				void describe(string configId) throws mdaq::utl::CmdException;
				void prepare(string configId) throws mdaq::utl::CmdException;
				void reset() throws mdaq::utl::CmdException;
				void configure(string configId, string dataRouter) throws mdaq::utl::CmdException;
				void start() throws mdaq::utl::CmdException;
				void stop() throws mdaq::utl::CmdException;
				void pause() throws mdaq::utl::CmdException;
				void resume() throws mdaq::utl::CmdException;
				void breakup() throws mdaq::utl::CmdException;
				void getStatus(out Status s) throws mdaq::utl::CmdException;
				idempotent string getDataLinks() throws mdaq::utl::CmdException;
			};
			
			interface ConfigurationManager
			{
				string getConfigIds() throws mdaq::utl::CmdException;
				string getConfigSet(string configId) throws mdaq::utl::CmdException;
				string getSubConfigIds(string type) throws mdaq::utl::CmdException;
				string addSubConfig(string type, string subConfig, string subConfigId) throws mdaq::utl::CmdException;
				string getSubConfig(string type, string subConfigId) throws mdaq::utl::CmdException;
			};
			
			interface PLLManager
			{
				void loadPLLConfig(float ckwFreq) throws mdaq::utl::CmdException;
				void writePLLRegister(long value) throws mdaq::utl::CmdException;
				void resetPLL() throws mdaq::utl::CmdException;
			};
			
			interface CoBoMutantLinkManager
			{
				void setupCoBoMutantLink(get::cobo::MutantLinkMode mode) throws mdaq::utl::CmdException;
			};
			
			interface AsAdPulserMgr
			{
				void setAsAdPulserDefaultVoltage(long voltage) throws mdaq::utl::CmdException;
				void triggerAsAdPulser(long voltage) throws mdaq::utl::CmdException;
				void configureAsAdExternalPulser(bool enabled, long voltage) throws mdaq::utl::CmdException;
				void configureAsAdPeriodicPulser(long voltage, long period) throws mdaq::utl::CmdException;
				void setRandomAsAdPulserEnabled(bool enabled) throws mdaq::utl::CmdException;
				void configureAsAdDoublePulserMode(long voltage) throws mdaq::utl::CmdException;
				void startAsAdPeriodicPulser() throws mdaq::utl::CmdException;
				idempotent void stopAsAdPeriodicPulser() throws mdaq::utl::CmdException;
			};
	}; /* module rc */
	
	interface GetEcc extends ::mdaq::Ecc, ::get::rc::StateMachine, rc::ConfigurationManager, rc::PLLManager, rc::CoBoMutantLinkManager, rc::AsAdPulserMgr
	{
		void describeNode(string nodeType, string instanceId, string configPath) throws mdaq::utl::CmdException;
		void prepareNode(string nodeType, string instanceId, string configPath) throws mdaq::utl::CmdException;
		void configureNode(string nodeType, string instanceId, string configPath) throws mdaq::utl::CmdException;

		void asadSetISPM(string moduleName, string signalNameISP1, string signalNameISP2)
			throws mdaq::utl::CmdException;
			
		void coboSwitchToSlowControlMode(bool slowMode) throws mdaq::utl::CmdException;
		
		void coboSelectSingleAgetReadoutChannel(long channelId) throws mdaq::utl::CmdException;
		
		void coboSelectSingleAgetTestChannel(long channelId) throws mdaq::utl::CmdException;
		
		void coboSetAgetChannelLsbThreshold(long channelId, long threshold) throws mdaq::utl::CmdException;
				
		void coboSetReadingClockFrequency(float frequency)
			throws mdaq::utl::CmdException;
			
		void coboSetWritingClockFrequency(float frequency, bool pllConfiEnabled, cobo::PllRegisterMap registers)
			throws mdaq::utl::CmdException;

		void setAsAdAlarmMonitoringEnabled(bool enabled) throws mdaq::utl::CmdException;
		void subscribeAlarms(bool enabled) throws mdaq::utl::CmdException;
		mt::AlarmService* alarmService() throws mdaq::utl::CmdException;
		
		cobo::LedManager* ledManager() throws mdaq::utl::CmdException;
		
		mutant::AlignmentServer* mutantAlignmentServer() throws mdaq::utl::CmdException;		
	};
}; /* module get */

#endif /* get_GetEcc_ice_INCLUDED */
