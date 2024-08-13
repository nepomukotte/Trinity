
#ifndef get_cobo_LedManager_ice_INCLUDED
#define get_cobo_LedManager_ice_INCLUDED

#include "mdaq/hw/Control.ice"

module get
{
	module cobo
	{
		/** Names of the 4 LEDS per AsAd (P, U and D are green, B is red). **/
		enum LedType { LedP, LedB, LedU, LedD, NumLedsPerAsAd };
		
		/** States available for each LED. **/
		enum LedState { LedOff, LedOn, LedSlowPulse, LedFastPulse };

		/** Interface for managing CoBo LED device **/
		interface LedManager
		{
			void blinkLEDs();
			void setLED(int led, bool enabled);
			void setLEDs(bool enabled);
			void flipLED(int led);
			void flipLEDs();
			void pulseLED(int led, int period);
			void modifyLED(LedType type, int asadId, LedState state);
		};
	};
};

#endif /* get_cobo_LedManager_ice_INCLUDED */