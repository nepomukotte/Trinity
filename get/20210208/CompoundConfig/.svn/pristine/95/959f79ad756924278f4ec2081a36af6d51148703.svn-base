// =============================================================================
// BinaryFormatFactory.h
// -----------------------------------------------------------------------------
// $Id$
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Patrick Sizun (patrick.sizun@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef CCfg_Io_BinaryFormatFactory_h_INCLUDED
#define CCfg_Io_BinaryFormatFactory_h_INCLUDED

#include "AbstractFormatFactory.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class creates BinaryFormat objects from URLs.
 *
 * @author sizun
 */
class BinaryFormatFactory : public AbstractFormatFactory
{
	static const int BINARY_EXT_COUNT = 3;
	static const char* const BINARY_EXTENSIONS[BINARY_EXT_COUNT];

public:
	~BinaryFormatFactory() {}
	std::auto_ptr<Format> createFormat(const Url& url) const;
};

}
}

#endif // CCfg_Io_BinaryFormatFactory_h_INCLUDED
