// =============================================================================
// IntervalGrammar.h
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

#ifndef CCfg_Grammar_IntervalGrammar_h_INCLUDED
#define CCfg_Grammar_IntervalGrammar_h_INCLUDED

#include "../Interval.h"
#include "ValueGrammar.h"
#include <boost/shared_ptr.hpp>
#include <boost/version.hpp>

namespace CCfg
{
namespace Grammar
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

#if BOOST_VERSION < 105300
typedef boost::shared_ptr< CCfg::Interval > IntervalPtr;
#else
struct IntervalPtr : public boost::shared_ptr< CCfg::Interval >
{
	typedef boost::detail::sp_dereference< CCfg::Interval >::type reference;
};
#endif

/**
 * This class is a closure (a parsing local context) for an interval parser.
 */
struct IntervalClosure : boost_spirit_classic::closure<IntervalClosure, IntervalPtr >
{
	member1 interval;
};

struct BoundModeClosure : boost_spirit_classic::closure<BoundModeClosure, Interval::BoundMode>
{
	member1 mode;
};

/**
 * This class defines the grammar for parsing the string representation of a CCfg::Interval.
 *
 * @author fchateau
 */
struct IntervalGrammar : boost_spirit_classic::grammar<IntervalGrammar, IntervalClosure::context_t>
{
	template <typename ScannerT>
	struct definition
	{
		definition(const Grammar::IntervalGrammar& self);
		const boost_spirit_classic::rule<ScannerT>& start() const;

		Grammar::ValueGrammar value_p;
		boost_spirit_classic::rule<ScannerT, BoundModeClosure::context_t> interval_begin;
		boost_spirit_classic::rule<ScannerT, BoundModeClosure::context_t> interval_end;
		boost_spirit_classic::rule<ScannerT> interval_p;
	};

	static void fromString(CCfg::Interval& interval, const std::string& str);
};

}
}

#endif // CCfg_Grammar_IntervalGrammar_h_INCLUDED
