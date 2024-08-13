///
///
/// @brief To change the char case in a string
///
/// @file change_case.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 22 feb 2022 
///

#ifndef CHANGE_CASE_HPP_
#define CHANGE_CASE_HPP_

#include <algorithm> // for std::transform
#include <cctype> // for std::tolower
#include <utility> // for std::forward
#include <string>
#include <string_view>

namespace fcutils{
namespace string{

struct to_lower_type{
	// Remember tolower, converts int. So to convert char you need to first cast it to unsigned char. An implicit onversion will tacke care of it
	auto operator()(unsigned char c){ return std::tolower(c); }
};

struct to_upper_type{
	// Remember tolower, converts int. So to convert char you need to first cast it to unsigned char. An implicit onversion will tacke care of it
	auto operator()(unsigned char c){ return std::toupper(c); }
};


auto to_lower(std::string s){std::transform( std::begin(s), std::end(s), std::begin(s), to_lower_type{}); return s;};
auto to_upper(std::string s){std::transform( std::begin(s), std::end(s), std::begin(s), to_upper_type{}); return s;};
auto to_lower(std::string_view c){std::string s{c}; return to_lower(std::move(s));};
auto to_upper(std::string_view c){std::string s{c}; return to_upper(std::move(s));};
auto to_lower(const char *c){std::string s{c}; return to_lower(std::move(s));};
auto to_upper(const char *c){std::string s{c}; return to_upper(std::move(s));};


} // namesapce string
} // namespae fcutils


#endif // CHANGE_CASE_HPP_ 
