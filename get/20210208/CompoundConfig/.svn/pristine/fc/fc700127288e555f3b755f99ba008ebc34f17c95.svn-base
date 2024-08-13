// =============================================================================
// EnumGrammar.hpp
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

#ifndef CCfg_Grammar_EnumGrammar_hpp_INCLUDED
#define CCfg_Grammar_EnumGrammar_hpp_INCLUDED

#include "EnumGrammar.h"
#include <complex>
#include <boost/spirit/include/phoenix1.hpp>

namespace CCfg
{
namespace Grammar
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Grammar::EnumGrammar::definition class methods body
////////////////////////////////////////////////////////////////////////////////

static void allocEnum(boost::shared_ptr<CCfg::Enum>& ptr);

void allocEnum(boost::shared_ptr<CCfg::Enum>& ptr)
{
	ptr = boost::shared_ptr<CCfg::Enum>(new CCfg::Enum);
}

/**
 * Constructs the grammar definition.
 */
template <typename ScannerT>
EnumGrammar::definition<ScannerT>::definition(const Grammar::EnumGrammar& self)
	: value_p(true, CCfg::Enum::RESERVED_CHARACTERS)
{
	using boost_spirit_classic::ch_p;
	using phoenix::arg1;
	using phoenix::bind;

	enum_p
		=	ch_p(CCfg::Enum::OPEN_CHAR) [bind(&allocEnum)(self.enumeration)]
		>>	!(value_p	[*self.enumeration = arg1]
			>> *(ch_p(CCfg::Enum::SEPARATOR) >> value_p	[*self.enumeration += arg1])
			)
		>> ch_p(CCfg::Enum::CLOSE_CHAR);
}

/**
 * Returns the grammar entry point.
 */
template <typename ScannerT>
const boost_spirit_classic::rule<ScannerT>&
EnumGrammar::definition<ScannerT>::start() const
{
	return enum_p;
}

}
}

#endif // CCfg_Grammar_EnumGrammar_hpp_INCLUDED
