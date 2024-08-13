
#ifndef get_cobo_CtrlNode_ice_INCLUDED
#define get_cobo_CtrlNode_ice_INCLUDED

#include "mdaq/hw/Control.ice"
#include "get/mt/AlarmService.ice"
#include "get/cobo/LedManager.ice"

module get
{
	module cobo
	{
	
		/** Interface for management of the AsAd voltage generator. **/
		interface AsAdPulserMgr
		{
			/** Injects charge in mV **/
			void setDefaultPulserVoltage(long voltage) throws mdaq::utl::CmdException;
			void triggerPulser(long voltage) throws mdaq::utl::CmdException;
			void configureExternalPulser(bool enabled, long voltage) throws mdaq::utl::CmdException;
			void configurePeriodicPulser(long voltage, long period) throws mdaq::utl::CmdException;
			void setRandomPulserEnabled(bool enabled) throws mdaq::utl::CmdException;
			void configureDoublePulserMode(long voltage) throws mdaq::utl::CmdException;
			void startPeriodicPulser() throws mdaq::utl::CmdException;
			idempotent void stopPeriodicPulser();
		};
		
		/** Interface for managing the monitoring of the AsAd alarms. **/
		interface AsAdAlarmMonitor
		{
			void setAsAdAlarmMonitoringEnabled(bool enabled) throws mdaq::utl::CmdException;
		};
		
		/** Interface for controlling a CoBo board **/
		interface CtrlNode extends ::mdaq::hw::Node, AsAdPulserMgr , AsAdAlarmMonitor, LedManager
		{
			mt::AlarmService* getAlarmService() throws mdaq::utl::CmdException;
		};

	};
};

#endif /* get_cobo_CtrlNode_ice_INCLUDED */