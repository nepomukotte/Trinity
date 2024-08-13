///
///
/// @brief A function to convert a char array into a string of hexadecimal char representation
///
/// @file char_to_hex.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 9 dic 2021 
///

#ifndef CHAR_TO_HEX_HPP_
#define CHAR_TO_HEX_HPP_
#include <iostream>
#include <string>
#include <charconv> // for std::to_chars
#include <system_error> // for std::errc and std::make_error_code
#include <cstdint> // for uint_fast16_t
#include <array>
#include <cstddef> // for std::byte
#include <cstring> // for std::memcpy
#include <type_traits> // for is_pointer

namespace fcutils{
namespace string{

/// The case with iterator. In this case is user responsibility to create the buffer and terminate after the int_to_hex call.
template <typename Char , typename Iter = char *> auto char_to_hex( Iter start, const Char * c, std::size_t length, const std::string & prefix="0x" ){
    std::array<Char, 10> temp;
    int base{16}; // The base, hexadecimal in this case
    // Copy the prefix
    for (auto & c: prefix) *start++=c;
    std::size_t diff{0};
    for(std::size_t i=0; i<length; ++i){
    	auto val=std::to_integer<uint_fast16_t>(static_cast<std::byte>(static_cast<unsigned char>( *(c+i) ) ) );
    	// Skip the leading zero, if the numerical value is less than base
    	std::size_t offset=val<base?1:0;
    	// set a leading zero, will be overwirtten in case the value is larger than base
    	temp[0]='0';
    	// Use if with inizialization
    	if(auto [ptr, ec] = std::to_chars(temp.data()+offset, temp.data() + temp.size(),val,base);
    	ec == std::errc()){
            diff=(ptr-temp.data());
            if constexpr (std::is_pointer<Iter>::value) {
            	std::memcpy(start,temp.data(),diff);
            	start+=diff;
            } else {
            	auto source=temp.data();
            	for(int i=0; i<diff; ++i) *start++=*source++;
            }
    	} else {
    		std::cout << " Error converting char: "<< std::make_error_code(ec).message() << std::endl;
    	}
    }
    // Return the first available iterator
    return start;
}

template <typename Char > std::string char_to_hex( const Char * c, std::size_t length, const std::string & prefix="0x"){
    std::array<Char, 10> temp;
    int base{16}; // The base, hexadecimal in this case
    std::string result{prefix};
    for(std::size_t i=0; i<length; ++i){
    	auto val=std::to_integer<uint_fast16_t>(static_cast<std::byte>(static_cast<unsigned char>( *(c+i) ) ) );
    	// Skip the leading zero, if the numerical value is less than base
    	std::size_t offset=val<base?1:0;
    	// set a leading zero, will be overwirtten in case the value is larger than base
    	temp[0]='0';
    	// Use if with inizialization
    	if(auto [ptr, ec] = std::to_chars(temp.data()+offset, temp.data() + temp.size(),val,base);
    	ec == std::errc()){
    		result.append(temp.data(),ptr);
    	} else {
    		std::cout << " Error converting char: "<< std::make_error_code(ec).message() << std::endl;
    	}
    }
    return result;
}

} // namespace string
} // namespace fcutils



#endif // CHAR_TO_HEX_HPP_ 
