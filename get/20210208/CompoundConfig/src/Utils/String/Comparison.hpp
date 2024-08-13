// =============================================================================
// Comparison.hpp
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

#ifndef Utils_String_Comparison_hpp_INCLUDED
#define Utils_String_Comparison_hpp_INCLUDED

#include <functional>
#include <locale>

namespace Utils
{
namespace String
{

/**
 * Case insensitive less functor for characters.
 */
class LessNoCaseChar
{
public:
	LessNoCaseChar(const std::ctype<char>& c) : ct(c) {}
	bool operator()(char x, char y) const
	{
		return ct.toupper(x) < ct.toupper(y);
	}
private:
	const std::ctype<char>& ct;
};

/**
 * Case insensitive equal functor for characters.
 */
class EqualNoCaseChar
{
public:
	EqualNoCaseChar(const std::ctype<char>& c) : ct(c) {}
	bool operator()(char x, char y) const
	{
		return ct.toupper(x) == ct.toupper(y);
	}
private:
	const std::ctype<char>& ct;
};

/**
 * Case insensitive less functor for strings.
 */
class LessNoCaseStr : public std::binary_function<std::string, std::string, bool>
{
public:
	LessNoCaseStr(const std::locale& L = std::locale::classic())
		: loc(L), ct(std::use_facet<std::ctype<char> >(loc))
	{
	}
	bool operator()(const std::string& x, const std::string& y) const
	{
		return std::lexicographical_compare(x.begin(), x.end(),
					y.begin(), y.end(),
					LessNoCaseChar(ct));
	}
private:
	std::locale loc;
	const std::ctype<char>& ct;
};

/**
 * Case insensitive equals functor for strings.
 */
class EqualNoCaseStr : public std::binary_function<std::string, std::string, bool>
{
public:
	EqualNoCaseStr(const std::locale& L = std::locale::classic())
		: loc(L), ct(std::use_facet<std::ctype<char> >(loc))
	{
	}
	bool operator()(const std::string& x, const std::string& y) const
	{
		if(x.size() != y.size())
		{
			return false;
		}
		std::string::const_iterator xCur = x.begin(), xEnd = x.end();
		std::string::const_iterator yCur = y.begin(), yEnd = y.end();
		EqualNoCaseChar cmp(ct);
		while(xCur != xEnd and yCur != yEnd)
		{
			if(not cmp(*xCur, *yCur))
			{
				return false;
			}
			++xCur;
			++yCur;
		}
		return true;
	}
private:
	std::locale loc;
	const std::ctype<char>& ct;
};

}
}


#endif // Utils_String_Comparison_hpp_INCLUDED
