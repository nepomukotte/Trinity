// =============================================================================
// EnumGrammar.cpp
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

#include "EnumGrammar.hpp"

namespace CCfg
{
namespace Grammar
{
using namespace boost_spirit_classic;

////////////////////////////////////////////////////////////////////////////////
// CCfg::Grammar::EnumGrammar class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Binds a CCfg::Enum with the result of parsing a string.
 */
void EnumGrammar::fromString(CCfg::Enum& enumeration, const std::string& str)
{
	boost::shared_ptr<CCfg::Enum> e;
	Grammar::EnumGrammar g;
	boost_spirit_classic::parse(str.c_str(), g[assign_a(e)], space_p);
	if(e != 0)
	{
		enumeration = *e;
	}
}

// explicitly instanciates the "definition" inner template class for standard scanner.
template struct EnumGrammar::definition<
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
