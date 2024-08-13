///
///
/// @brief To convert an ASCII representation in hexadecimal values into a char
///
/// @file hex_to_char.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 6 feb 2022 
///

#ifndef HEX_TO_CHAR_HPP_
#define HEX_TO_CHAR_HPP_
#include <iostream>
#include <string>
#include <charconv> // for std::from_chars
#include <system_error> // for std::errc and std::make_error_code
#include <cstdint> // for uint_fast16_t
#include <tuple>
#include <cstddef> // for std::byte
#include <cstring> // for std::memcpy
#include <type_traits> // for is_pointer


namespace fcutils{
namespace string{

template <typename Char, typename Iter=char * > auto hex_to_char( Iter t, const Char * c, std::size_t length){
    int base{16}; // The base, hexadecimal in this case
    const std::size_t BYTE_SIZE{2};
    // Skip the hexadecimal prefix if any
    auto start=c;
    if( (c[0]=='0') && ((c[1]=='x') || (c[1]=='X')) ) start+=BYTE_SIZE;
    uint_fast8_t val{};
    while(start<(c+length)){
    	if(auto [ptr, ec] = std::from_chars(start, start+BYTE_SIZE,val,base);
    	    	ec == std::errc()){
    				start=ptr;
    				// leveraging relaxed inizialization in C++17
    				std::byte temp{val};
    	    		if constexpr ( std::is_same< typename std::remove_reference<Iter>::type, std::string>::value ) {
    	    			t+=static_cast<unsigned char>(temp);
    	    		} else {
    	    			*t++=static_cast<unsigned char>(temp);
    	    		}
    	    	} else {
    	    		std::cout << " Error converting from hex: "<< std::make_error_code(ec).message() << std::endl;
    	    	}
    } // while
    return std::tuple<Iter, const Char *>{t,start};
}

template <typename Char > std::string hex_to_char( const Char * c, std::size_t length){
    int base{16}; // The base, hexadecimal in this case
    std::string result;
    result.reserve(length);
    auto [ptr,str] = hex_to_char<Char, std::string &>( result, c,length);
    return result;
}

} // namespace string
} // namespace fcutils



#endif // HEX_TO_CHAR_HPP_ 
