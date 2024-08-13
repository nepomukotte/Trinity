// =============================================================================
// ValueGrammar.h
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

#ifndef CCfg_Grammar_ValueGrammar_h_INCLUDED
#define CCfg_Grammar_ValueGrammar_h_INCLUDED

#include "../Value.h"
#include <boost/version.hpp>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_closure.hpp>
namespace boost_spirit_classic = boost::spirit::classic;

namespace CCfg
{
namespace Grammar
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class is a closure (a parsing local context) for a CCfg::ValueGrammar parser.
 */
struct ValueClosure : boost_spirit_classic::closure<ValueClosure, Value>
{
    member1 val;
};

/**
 * This class is a closure (a parsing local context) for a string parser.
 */
struct StringClosure : boost_spirit_classic::closure<StringClosure, std::string>
{
	member1 str;
};

/**
 * This class is a closure (a parsing local context) for a boolean parser.
 */
struct BoolClosure : boost_spirit_classic::closure<BoolClosure, bool>
{
	member1 boolean;
};

/**
 * This class defines the grammar for parsing the string representation of a CCfg::ValueGrammar.
 *
 * @author fchateau
 */
struct ValueGrammar : boost_spirit_classic::grammar<ValueGrammar, ValueClosure::context_t>
{
	template <typename ScannerT>
	struct definition
	{
		definition(const Grammar::ValueGrammar& self);
		const boost_spirit_classic::rule<ScannerT>& start() const;

		boost_spirit_classic::rule<ScannerT> test_separator;
		boost_spirit_classic::rule<ScannerT, BoolClosure::context_t> bool_p;
		boost_spirit_classic::rule<ScannerT, StringClosure::context_t> escape_str_p;
		boost_spirit_classic::rule<ScannerT, StringClosure::context_t> string_p;
		boost_spirit_classic::rule<ScannerT> value_p;
	};

	ValueGrammar(bool extractQuoted, const char* restrictedChars);

	static void fromString(CCfg::Value& val, const std::string& str, const char* restrictedChars = "");
	static void popTail(std::string& str);

	bool extractQuoted;
	const char* restrictedChars;
};

}
}

#endif // CCfg_Grammar_ValueGrammar_h_INCLUDED
