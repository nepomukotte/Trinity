// =============================================================================
// RangeGrammar.cpp
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

#include "RangeGrammar.hpp"

namespace CCfg
{
namespace Grammar
{
using namespace boost_spirit_classic;

////////////////////////////////////////////////////////////////////////////////
// CCfg::Grammar::RangeGrammar class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Binds a CCfg::Enum with the result of parsing a string.
 */
std::auto_ptr<CCfg::Range> RangeGrammar::fromString(const std::string& str)
{
	boost::shared_ptr<CCfg::Range> r;
	Grammar::RangeGrammar g;
	boost_spirit_classic::parse(str.c_str(), g[assign_a(r)], space_p);

	return r != 0 ? r->clone() : std::auto_ptr<CCfg::Range>();
}

// explicitly instanciates the "definition" inner template class for standard scanner.
template struct RangeGrammar::definition<
	scanner<
		char const*,
		scanner_policies<
			skipper_iteration_policy<iteration_policy>,
			match_policy,
			action_policy
		>
	>
>;

}
}
