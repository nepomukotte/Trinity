///
///
/// @brief A new file
///
/// @file int_to_hex.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 23 gen 2022 
///

#ifndef INT_TO_HEX_HPP_
#define INT_TO_HEX_HPP_
#include <iostream>
#include <string>
#include <charconv> // for std::to_chars
#include <system_error> // for std::errc and std::make_error_code
#include <cstring> // for std::memcpy
#include <array>
#include <type_traits> // for is_pointer

namespace fcutils{
namespace string{

/// The case with iterator. In this case is user responsibility to create the buffer and terminate after the int_to_hex call.
template <typename Int , typename Iter = char *> auto int_to_hex( Iter start, const Int & val, const std::string & prefix="0x", bool pad=true ){
    constexpr const std::size_t n_char=sizeof(Int)*2;
    std::array<char, n_char> temp;
    int base{16}; // The base, hexadecimal in this case
    std::size_t diff{0};
    if(auto [ptr, ec] = std::to_chars(temp.data(), temp.data() + n_char ,val,base); ec == std::errc()){
        diff=(ptr-temp.data());
        int n_pad=n_char-diff;
        // Copy the prefix
        for (auto & c: prefix) *start++=c;
        // Fill with zeros if padding is requested
        if(pad && (n_pad>0) ) for(int i=0; i<n_pad; ++i) *start++='0';
        if constexpr (std::is_same<Iter, char *>::value || std::is_same<Iter, unsigned char *>::value) {
        	// If is a pointer to char, just use memcpy
        	std::memcpy(start,temp.data(),diff);
        	start+=diff;
        } else {
        	// Is an iterator, just postincrement and dereference.
        	auto source=temp.data();
        	for(int i=0; i<diff; ++i) *start++=*source++;
        }
    } else {
        std::cout << " Error converting: "<< std::make_error_code(ec).message() << std::endl;
    }
    // return the first available pointer
    return start;
}

template <typename Int > std::string & int_to_hex( std::string & start, const Int & val, const std::string & prefix="0x", bool pad=true ){
    constexpr const std::size_t n_char=sizeof(Int)*2;
    std::array<char, n_char> temp;
    int base{16}; // The base, hexadecimal in this case
    std::size_t diff{0};
    if(auto [ptr, ec] = std::to_chars(temp.data(), temp.data() + n_char ,val,base); ec == std::errc()){
        diff=(ptr-temp.data());
        int n_pad=n_char-diff;
        // If is a string, then append
        // Copy the prefix
        start+=prefix;
        // Fill with zeros if padding is requested
        if(pad && (n_pad>0) ) start.append(n_pad,'0');
        auto source=temp.data();
        start.append(temp.data(),ptr);
    } else {
        std::cout << " Error converting: "<< std::make_error_code(ec).message() << std::endl;
    }
    // return the first available pointer
    return start;
}

template <typename Int > std::string int_to_hex( const Int & val, const std::string & prefix="0x", bool pad=true ){
    constexpr const std::size_t n_char=sizeof(Int)*2;
    std::string result;
    result.reserve(n_char+prefix.size());
    return int_to_hex(result,val, prefix, pad);
}
//template <typename Int > std::string int_to_hex( const Int & val, const std::string & prefix="0x", bool pad=true ){
//    constexpr const std::size_t n_char=sizeof(Int)*2;
//    std::array<char, n_char> temp;
//    int base{16}; // The base, hexadecimal in this case
//    std::string result;
//    result.reserve(n_char+prefix.size());
//    result=prefix;
//    if(auto [ptr, ec] = std::to_chars(temp.data(), temp.data() + temp.size(),val,base); ec == std::errc()){
//        std::size_t diff=(ptr-temp.data());
//        int n_pad=n_char-diff;
//        // Fill with zeros if padding is requested
//        if(pad && (n_pad>0) ) result.append(n_pad,'0');
//        result.append(temp.data(),ptr);
//    } else {
//        std::cout << " Error converting: "<< std::make_error_code(ec).message() << std::endl;
//    }
//    return result;
//}


} // namespace string
} // namesapce fcutils

#endif // INT_TO_HEX_HPP_
