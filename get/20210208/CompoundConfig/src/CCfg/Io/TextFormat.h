// =============================================================================
// TextFormat.h
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

#ifndef CCfg_Io_TextFormat_h_INCLUDED
#define CCfg_Io_TextFormat_h_INCLUDED

#include "Format.h"

namespace CCfg
{
namespace Io
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * An intermediary interface to distinguish format families.
 *
 * @author fchateau
 */
class TextFormat : public Format
{
public:
	TextFormat();
	virtual ~TextFormat();

	virtual std::string getConfigName(Storage& storage);
};

}

}

#endif // CCfg_Io_TextFormat_h_INCLUDED
