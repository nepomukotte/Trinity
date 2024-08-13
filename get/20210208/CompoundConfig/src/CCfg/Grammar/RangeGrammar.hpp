// =============================================================================
// RangeGrammar.hpp
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

#ifndef CCfg_Grammar_RangeGrammar_hpp_INCLUDED
#define CCfg_Grammar_RangeGrammar_hpp_INCLUDED

#include "RangeGrammar.h"
#include <complex>
#include <boost/spirit/include/phoenix1.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>


namespace CCfg
{
namespace Grammar
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Grammar::RangeGrammar::definition class methods body
////////////////////////////////////////////////////////////////////////////////


/**
 * Constructs the grammar definition.
 */
template <typename ScannerT>
RangeGrammar::definition<ScannerT>::definition(const Grammar::RangeGrammar& self)
{
	using phoenix::arg1;
	using phoenix::bind;

	range_p
		= enum_p		[self.range_ptr = arg1]
		| interval_p	[self.range_ptr = arg1];
}

/**
 * Returns the grammar entry point.
 */
template <typename ScannerT>
const boost_spirit_classic::rule<ScannerT>&
RangeGrammar::definition<ScannerT>::start() const
{
	return range_p;
}

}
}

#endif // CCfg_Grammar_RangeGrammar_hpp_INCLUDED
