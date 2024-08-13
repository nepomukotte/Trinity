// =============================================================================
// EnumGrammar.h
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

#ifndef CCfg_Grammar_EnumGrammar_h_INCLUDED
#define CCfg_Grammar_EnumGrammar_h_INCLUDED

#include "../Enum.h"
#include "ValueGrammar.h"
#include "boost/shared_ptr.hpp"
#include <boost/version.hpp>

namespace CCfg
{
namespace Grammar
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

#if BOOST_VERSION < 105300
typedef boost::shared_ptr< CCfg::Enum > EnumPtr;
#else
struct EnumPtr : public boost::shared_ptr< CCfg::Enum >
{
	typedef boost::detail::sp_dereference< CCfg::Enum >::type reference;
};
#endif

/**
 * This class is a closure (a parsing local context) for an enumeration parser.
 */
struct EnumClosure : boost_spirit_classic::closure<EnumClosure, EnumPtr >
{
	member1 enumeration;
};

/**
 * This class defines the grammar for parsing the string representation of a CCfg::EnumGrammar.
 *
 * @author fchateau
 */
struct EnumGrammar : boost_spirit_classic::grammar<EnumGrammar, EnumClosure::context_t>
{
	template <typename ScannerT>
	struct definition
	{
		definition(const Grammar::EnumGrammar& self);
		const boost_spirit_classic::rule<ScannerT>& start() const;

		Grammar::ValueGrammar value_p;
		boost_spirit_classic::rule<ScannerT> enum_p;
	};

	static void fromString(CCfg::Enum& enumeration, const std::string& str);
};

}
}


#endif // CCfg_Grammar_EnumGrammar_h_INCLUDED
