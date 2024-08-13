///
///
/// @brief An anonymous namespace to define i/o stream operator for std::byte
///
/// @file byte_iostream.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 27 feb 2022 
///

#ifndef UTILITY_BYTE_IOSTREAM_HPP_
#define UTILITY_BYTE_IOSTREAM_HPP_
#include <iosfwd>
//#include <bitset>
#include <cstddef> // for stb::byte
#include <limits>

namespace {

std::istream& operator>> (std::istream & i,  std::byte & b)
{
//	// read into a bitset:
//	std::bitset<std::numeric_limits<unsigned char>::digits> bs;
//	i >> bs;
//	// without failure, convert to std::byte:
//	if (! i.fail()) {
//		b = static_cast<std::byte>(bs.to_ulong()); // OK
//	}
	unsigned int ui;
	i >> ui;
	if (! i.fail()) {
		b = static_cast<std::byte>(ui); // OK
	}
	return i;
}

std::ostream & operator<<(std::ostream & o, const std::byte b ){
	return o << std::to_integer<unsigned int>(b) ;
}

}




#endif // UTILITY_BYTE_IOSTREAM_HPP_ 
