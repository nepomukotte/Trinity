// =============================================================================
// Config.h
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef CCfg_Io_Config_h_INCLUDED
#define CCfg_Io_Config_h_INCLUDED


#include "CCfg/Object.h"

namespace CCfg
{

struct IoTestConfig
{
	CCfg::Object daq;
	CCfg::Object config;
	CCfg::Object ambiguityTest;

	IoTestConfig();
	~IoTestConfig();
};

}

#endif // CCfg_Io_Config_h_INCLUDED
