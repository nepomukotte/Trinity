// $Id$
#ifndef get_Mutant_ice_INCLUDED
#define get_Mutant_ice_INCLUDED

module get
{
	module mutant
	{
	
		/** Clock network **/
		enum Clock { GMC, WSCA };
		
		/** Results of alignment procedure for a single slot **/
		struct SlotAlignmentResult
		{
			/** True if AMC present **/
			bool isAmcPresent;
			/** Delay value in ps **/
			int amcDeltaDelay;
			/** Value programmed in Delay_line register **/
			int amcProgDelay;
		};
		
		/** Results of alignment procedure, mapped by AMC slot identifier **/
		dictionary<short, SlotAlignmentResult> AlignmentResults;
		
		/** Alarm thrown by alignment procedure. **/
		exception AlignmentAlarm
		{
			string reason;
		};
		
		/** Results of second step of multi-shelve alignment procedure**/
		struct MultiShelveAlignmentResults
		{
			int measuredDelay;
			int programmedDelay;
			string comment;
		};
		
		/** Interface for Mutant alignment procedure **/
		interface AlignmentServer
		{
			void alignDelayTDC(Clock c, bool verbose, out AlignmentResults results)
				throws AlignmentAlarm;
			void multiShelveAlignDelayStep1(Clock c, bool verbose, out MultiShelveAlignmentResults results)
				throws AlignmentAlarm;
			void multiShelveAlignDelayStep2(Clock c, bool verbose, out MultiShelveAlignmentResults results)
				throws AlignmentAlarm;
			void multiShelveAlignDelayStep3(Clock c, bool verbose, out MultiShelveAlignmentResults results)
				throws AlignmentAlarm;
		};
	};
};

#endif /* get_Mutant_ice_INCLUDED */
