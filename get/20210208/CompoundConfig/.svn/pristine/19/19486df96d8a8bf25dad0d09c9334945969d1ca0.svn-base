// =============================================================================
// FormatFactory.cpp
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

#include "FormatFactory.h"
#if BUILD_CCFG_FMT_XML
#include "XmlFormatFactory.h"
#endif
#if BUILD_CCFG_FMT_BINARY
#include "BinaryFormatFactory.h"
#endif
#include "UnhandledUrl.h"
#include "Format.h"
#include "Url.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::FormatFactory class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Private constructor.
 * Automatically registers the Xml format and Binary format factories if available.
 */
FormatFactory::FormatFactory()
{
#if BUILD_CCFG_FMT_XML
	registerFormat(std::auto_ptr<AbstractFormatFactory>(new XmlFormatFactory));
#endif
#if BUILD_CCFG_FMT_BINARY
	registerFormat(std::auto_ptr<AbstractFormatFactory>(new BinaryFormatFactory));
#endif
}

/*
 * (see inherited doc)
 */
std::auto_ptr<Format> FormatFactory::createFormat(const Url& url) const
{
	std::auto_ptr<Format> format;
	for(boost::ptr_vector<AbstractFormatFactory>::const_iterator curFactory = factories.begin(),
		endFactory = factories.end();
		curFactory != endFactory and format.get() == 0;
		++curFactory)
	{
		std::auto_ptr<Format> format2 = curFactory->createFormat(url);
		format = format2; // hack for vxWorks
	}
	if(format.get() == 0)
	{
		throw UnhandledUrl(url.format(), "Format");
	}
	return format;
}

/**
 * Registers a format factory to enable the recognition of a format in the whole
 * application.
 */
void FormatFactory::registerFormat(std::auto_ptr<AbstractFormatFactory> factory)
{
	factories.push_back(factory.release());
}

/**
 * Returns the singleton instance.
 */
FormatFactory& FormatFactory::instance()
{
	static FormatFactory inst;
	return inst;
}

}
}
