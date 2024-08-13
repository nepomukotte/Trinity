// =============================================================================
// RangeGrammar.h
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

#ifndef CCfg_Grammar_RangeGrammar_h_INCLUDED
#define CCfg_Grammar_RangeGrammar_h_INCLUDED

#include "../Range.h"
#include "EnumGrammar.h"
#include "IntervalGrammar.h"
#include <boost/shared_ptr.hpp>

namespace CCfg
{
namespace Grammar
{

/**
 * This class is a closure (a parsing local context) for a range parser.
 */
struct RangeClosure : boost_spirit_classic::closure<RangeClosure, boost::shared_ptr<CCfg::Range> >
{
	member1 range_ptr;
};

/**
 * This class defines the grammar for parsing the string representation of a CCfg::Enum.
 *
 * @author fchateau
 */
struct RangeGrammar : boost_spirit_classic::grammar<RangeGrammar, RangeClosure::context_t>
{
	template <typename ScannerT>
	struct definition
	{
		definition(const Grammar::RangeGrammar& self);
		const boost_spirit_classic::rule<ScannerT>& start() const;

		Grammar::EnumGrammar enum_p;
		Grammar::IntervalGrammar interval_p;
		boost_spirit_classic::rule<ScannerT> range_p;
	};

	static std::auto_ptr<CCfg::Range> fromString(const std::string& str);
};

}
}

#endif // CCfg_Grammar_RangeGrammar_h_INCLUDED

