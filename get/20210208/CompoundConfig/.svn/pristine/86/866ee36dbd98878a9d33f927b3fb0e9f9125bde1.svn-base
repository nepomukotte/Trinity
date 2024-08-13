// =============================================================================
// ValueGrammar.hpp
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

#ifndef CCfg_Grammar_ValueGrammar_hpp_INCLUDED
#define CCfg_Grammar_ValueGrammar_hpp_INCLUDED

#include "ValueGrammar.h"
#include <boost/spirit/include/classic_confix.hpp>
#include <boost/spirit/include/classic_escape_char.hpp>
#include <boost/spirit/include/classic_chset.hpp>
#include <boost/spirit/include/classic_if.hpp>
#include <complex>
#include <boost/spirit/include/phoenix1.hpp>

namespace CCfg
{
namespace Grammar
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Grammar::ValueGrammar::definition class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructs the grammar definition.
 */
template <typename ScannerT>
ValueGrammar::definition<ScannerT>::definition(const Grammar::ValueGrammar& self)
{
	using boost_spirit_classic::strict_real_p;
	using boost_spirit_classic::hex_p;
	using boost_spirit_classic::int_p;
	using boost_spirit_classic::str_p;
	using boost_spirit_classic::end_p;
	using boost_spirit_classic::if_p;
	using boost_spirit_classic::as_lower_d;
	using boost_spirit_classic::c_escape_ch_p;
	using boost_spirit_classic::chset_p;
	using boost_spirit_classic::nothing_p;
	using boost_spirit_classic::confix_p;
	using boost_spirit_classic::lexeme_d;
	using boost_spirit_classic::epsilon_p;
	using phoenix::arg1;
	using phoenix::bind;

	value_p
		= strict_real_p			[self.val = arg1]	>> test_separator	// ensures the number is followed by a valid separator.
		| lexeme_d["0x" >> hex_p[self.val = arg1]]	>> test_separator
		| int_p					[self.val = arg1]	>> test_separator
		| bool_p				[self.val = arg1]	>> test_separator
		| escape_str_p			[self.val = arg1]
		| string_p				[self.val = arg1]
		;
	bool_p
		= as_lower_d
		[
			str_p("true")		[bool_p.boolean = true]
			| str_p("false")	[bool_p.boolean = false]
		];
	escape_str_p
		= if_p(phoenix::bind(&ValueGrammar::extractQuoted)(var(self)))
		[
		 	confix_p('"',
		 		(
					*c_escape_ch_p	[escape_str_p.str += arg1] // concatenate the closure with each character
				)
				[phoenix::bind(ValueGrammar::popTail)(escape_str_p.str)],
				'"')
		]
		.else_p
		[
		 	nothing_p
		];
	string_p
		= lexeme_d
		[
			+(~chset_p(self.restrictedChars))	[string_p.str += arg1]
		];
	test_separator
		= epsilon_p(chset_p(self.restrictedChars) | end_p);
}

/**
 * Returns the grammar entry point.
 */
template <typename ScannerT>
const boost_spirit_classic::rule<ScannerT>&
ValueGrammar::definition<ScannerT>::start() const
{
	return value_p;
}

}
}

#endif // CCfg_Grammar_ValueGrammar_hpp_INCLUDED
