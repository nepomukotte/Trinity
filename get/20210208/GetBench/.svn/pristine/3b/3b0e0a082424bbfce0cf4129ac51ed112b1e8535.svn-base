// $Id$
#ifndef get_mt_AlarmService_ice_INCLUDED
#define get_mt_AlarmService_ice_INCLUDED


module get
{
	module mt
	{
		/** Callback interface for clients to handle alarms. **/
		interface AlarmCallback
		{
			void sendAsadAlarm(string coboId, long alarmBits);
		};
		
		interface AlarmService
		{
			/** Registers a subscriber **/
			idempotent void subscribe(string id, AlarmCallback* cb);
			/** Sends AsAd alarm o all subscribers **/
			void sendAsadAlarm(string coboId, long alarmBits);
			/** Unregisters a subscriber **/
			idempotent void unsubscribe(string id);
			/** Deletes all subscribers **/
			idempotent void reset();
		};
	};
};

#endif /* get_mt_AlarmService_ice_INCLUDED */