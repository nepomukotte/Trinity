// =============================================================================
// IntervalGrammar.hpp
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

#ifndef CCfg_Grammar_IntervalGrammar_hpp_INCLUDED
#define CCfg_Grammar_IntervalGrammar_hpp_INCLUDED

#include "IntervalGrammar.h"
#include <complex>
#include <boost/spirit/include/phoenix1.hpp>

namespace CCfg
{
namespace Grammar
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Grammar::IntervalGrammar::definition class methods body
////////////////////////////////////////////////////////////////////////////////

static void allocInterval(boost::shared_ptr<CCfg::Interval>& ptr);

void allocInterval(boost::shared_ptr<CCfg::Interval>& ptr)
{
	ptr = boost::shared_ptr<CCfg::Interval>(new CCfg::Interval);
}

/**
 * Constructs the grammar definition.
 */
template <typename ScannerT>
IntervalGrammar::definition<ScannerT>::definition(const Grammar::IntervalGrammar& self)
	: value_p(true, CCfg::Interval::RESERVED_CHARACTERS)
{
	using boost_spirit_classic::ch_p;
	using boost_spirit_classic::str_p;
	using phoenix::arg1;
	using phoenix::bind;

	interval_begin
		= ch_p(CCfg::Interval::LBOUND_INC)	[interval_begin.mode = CCfg::Interval::IncludeBound]
		| ch_p(CCfg::Interval::LBOUND_EXC)	[interval_begin.mode = CCfg::Interval::ExcludeBound]
		;
	interval_end
		= ch_p(CCfg::Interval::UBOUND_INC)	[interval_end.mode = CCfg::Interval::IncludeBound]
		| ch_p(CCfg::Interval::UBOUND_EXC)	[interval_end.mode = CCfg::Interval::ExcludeBound]
		;
	interval_p
		= interval_begin	[bind(&allocInterval)(self.interval)]
		                	[bind(&CCfg::Interval::setLowerBoundMode)(*self.interval, arg1)]
		>>	(	str_p(CCfg::Interval::MINUS_INFINITY)	[bind(&CCfg::Interval::setLowerBoundMode)(*self.interval, CCfg::Interval::Unbounded)]
			|	value_p		[bind(&CCfg::Interval::setLowerBound)(*self.interval, arg1)]
			)
		>> ch_p(CCfg::Interval::SEPARATOR)
		>> ((str_p(CCfg::Interval::PLUS_INFINITY)	[bind(&CCfg::Interval::setUpperBoundMode)(*self.interval, CCfg::Interval::Unbounded)]
			>>	interval_end)
			|
			(	value_p			[bind(&CCfg::Interval::setUpperBound)(*self.interval, arg1)]
			>>	interval_end	[bind(&CCfg::Interval::setUpperBoundMode)(*self.interval, arg1)]
			))
		;
}

/**
 * Returns the grammar entry point.
 */
template <typename ScannerT>
const boost_spirit_classic::rule<ScannerT>&
IntervalGrammar::definition<ScannerT>::start() const
{
	return interval_p;
}

}
}

#endif // CCfg_Grammar_IntervalGrammar_hpp_INCLUDED
