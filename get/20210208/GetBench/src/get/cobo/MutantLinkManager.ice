
#ifndef get_cobo_MutantLinkManager_ice_INCLUDED
#define get_cobo_MutantLinkManager_ice_INCLUDED

#include "mdaq/Ecc.ice"

module get
{
	module cobo
	{
		/** CoBo-Mutant link modes:
		  - Normal Operation (Input buffer enabled)
		  - CoBo Disconnected from Mutant network
		  - GMC Feedback Mode (Alignment)
		  - WSCA Feedback Mode (Alignment)
		**/
		enum MutantLinkMode { Normal, Disconnected, GMC, WSCA };
		
		/** Interface for setting the CoBo-Mutant link mode **/
		interface MutantLinkManager
		{
			void setupMutantLink(MutantLinkMode mode) throws mdaq::utl::CmdException;
		};

	};
};


#endif /* get_cobo_MutantLinkManager_ice_INCLUDED */