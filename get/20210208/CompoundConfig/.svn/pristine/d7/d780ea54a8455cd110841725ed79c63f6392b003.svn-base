// =============================================================================
// ValueGrammar.cpp
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

#include "ValueGrammar.hpp"

namespace CCfg
{
namespace Grammar
{
using namespace boost_spirit_classic;

////////////////////////////////////////////////////////////////////////////////
// CCfg::Grammar::ValueGrammar class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param extractQuoted indicates whether the parser should extract strings inside
 * quotes or keep the quotes.
 * @param restrictedChars the list of characters that are not
 * allowed in non-delimited strings. By default all characters are allowed,
 * but if the string begins by a quote, it will be considered as a delimited string.
 */
ValueGrammar::ValueGrammar(bool extractQuoted, const char* restrictedChars)
	: extractQuoted(extractQuoted), restrictedChars(restrictedChars)
{

}

/**
 * Binds a CCfg::ValueGrammar with the result of parsing a string.
 */
void ValueGrammar::fromString(CCfg::Value& val, const std::string& str, const char* restrictedChars)
{
	Grammar::ValueGrammar g(false, restrictedChars);
	boost_spirit_classic::parse(str.c_str(), g[assign_a(val)], space_p);
}

/**
 * Removes the last character of a string.
 * This is needed because the semantic action on "c_escape_ch_p" is
 * called on the last quote of a "escape_str_p". This is a limitation of
 * Spirit that may be fixed in a later version.
 */
void
ValueGrammar::popTail(std::string& str)
{
	if(str.size() > 0)
	{
		str.resize(str.size() - 1);
	}
}

// explicitly instanciates the "definition" inner template class for standard scanner.
template struct ValueGrammar::definition<
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
