///
///
/// @brief A new file
///
/// @file string_to_object.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 4 apr 2022 
///

#ifndef STRING_TO_OBJECT_HPP_
#define STRING_TO_OBJECT_HPP_
#include <string_view>
#include <iostream>

#include "fcutils/string/change_case.hpp"
#include "fcutils/string/hash.hpp"

namespace {
	using string_type=std::string_view;
	const std::string_view DEFAULTKEY("def");
	const std::size_t H_DEFAULTKEY{fcutils::string::hash(DEFAULTKEY.data())};
}


namespace fcutils{
namespace monitor{


template<bool Hash=true, class TranslationMap>  auto string_to_object(std::string key, const TranslationMap &m){
	key=fcutils::string::to_lower(key);
	std::size_t hkey{0ul};
	auto tt = m.end();
	auto last = m.end();
	if constexpr(Hash){
		// Use the string hash as map key
		hkey = fcutils::string::hash(key.data());
		tt=m.find(hkey);
	} else {
		// Use the std::string_view instead of the string hash
		tt=m.find(key);
	}
	if(tt==last){
		std::cout <<" Ooppsss! Key: " << key << " not found. "
				<<" Using the default value: " << DEFAULTKEY << std::endl;
		if constexpr(Hash){
			tt=m.find(hkey);
		} else {
			tt=m.find(key);
		}
		if(tt==last) {
			std::cout <<" Ooppsss! Default Key: " << DEFAULTKEY << " not found. Using an empty object." << std::endl;
		}
	}

	return tt->second;
}


} // namespace monitor
} // namespace fcutils




#endif // STRING_TO_OBJECT_HPP_ 
