///
///
/// @brief A constexpr hash algorithm-.
/// Thanks to Nicolai M. Josuttis: Programming with C++17, for the idea.
/// See the djb2 algorithm at http://www.cse.yorku.ca/~oz/hash.html for the source of this algorithm.
///
/// @file hash.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 29 mar 2022 
///

#ifndef HASH_HPP_
#define HASH_HPP_
#include <cstddef>

namespace fcutils{
namespace string{

// Since C++17 we can have constexpr lambda!
// Note that the first constepr refers to the closure type, i.e. hashed, the second, refers to the operator(). The latest is the one that matters for the lambda.
// Anyway, lambda operator are by default constexpr since C++17
	constexpr auto hash_xor_op = [](auto h, char c) constexpr {return h*33 ^ c;};
	constexpr auto hash_sum_op = [](auto h, char c) constexpr {return h*33 + c;};
	constexpr auto generic_hash = [](const char* str, auto combine ) constexpr {
		std::size_t hh = 5381; // initial hash value
		while (*str != '\0') {
			hh = combine(hh, *str++); // combine hash value with next character
		}
		return hh;
	};
	constexpr auto hash = [](const char* str) constexpr {
		std::size_t hh = 5381; // initial hash value
		while (*str != '\0') {
			hh = hh * 33 ^ *str++;
		}
		return hh;
	};

} // namesapce string
} // namesapce fcutils




#endif // HASH_HPP_ 
