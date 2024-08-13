// =============================================================================
// FormatFactory.h
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

#ifndef CCfg_Io_FormatFactory_h_INCLUDED
#define CCfg_Io_FormatFactory_h_INCLUDED

#include "AbstractFormatFactory.h"
#include <boost/ptr_container/ptr_vector.hpp>

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class is a singleton that creates formats from URLs using a list of
 * registered factories.
 * The Xml format is automatically registered.
 *
 * @author fchateau
 */
class FormatFactory : public AbstractFormatFactory
{
	boost::ptr_vector<AbstractFormatFactory> factories;

public:
	~FormatFactory() {}
	std::auto_ptr<Format> createFormat(const Url& url) const;
	void registerFormat(std::auto_ptr<AbstractFormatFactory> factory);

	static FormatFactory& instance();

private:
	FormatFactory();
	FormatFactory(const FormatFactory& rvalue);
	FormatFactory& operator=(const FormatFactory& rvalue);
};

}
}

#endif // CCfg_Io_FormatFactory_h_INCLUDED
