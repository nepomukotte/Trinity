///
///
/// @brief A new file
///
/// @file hex_to_int.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 7 feb 2022 
///

#ifndef HEX_TO_INT_HPP_
#define HEX_TO_INT_HPP_

#include <iostream>
#include <string>
#include <charconv> // for std::from_chars
#include <system_error> // for std::errc and std::make_error_code
#include <cstdint> // for uint_fast16_t
#include <tuple>
#include <cstddef> // for std::byte
#include <cstring> // for std::memcpy
#include <type_traits> // for is_pointer
#include <iterator> // for iterator_traits


namespace fcutils{
namespace string{

template< typename T, std::size_t ByteSize=1 > struct type_size{
	enum: std::size_t { value=sizeof(T)*ByteSize };
	using type = T;
	constexpr explicit operator std::size_t() const { return value; }
};

template <typename Int, std::size_t IntSize=sizeof(Int), bool DoCheck=true > auto hex_to_int( const char * c, Int & val){
	if(!c)  return std::tuple<const char *, Int>{nullptr,val};
    int base{16}; // The base, hexadecimal in this case
    const std::size_t BYTE_SIZE{2},INT_SIZE_IN_CHAR{2*IntSize};
    // Skip the hexadecimal prefix if any
    auto start=c;
    if constexpr ( DoCheck ) {
    	if( (c[0]=='0') && ((c[1]=='x') || (c[1]=='X')) ) start+=BYTE_SIZE;
    }
    if(auto [ptr, ec] = std::from_chars(start, start+INT_SIZE_IN_CHAR,val,base);
    ec == std::errc()){
    	start=ptr;
    } else {
    	std::cout << " Error converting from hex: "<< std::make_error_code(ec).message() << std::endl;
    	start=nullptr;
    }
    return std::tuple<const char *, Int>{start,val};
}

template <typename Iter, typename SizeType = std::size_t, bool DoCheck=true > auto hex_to_int( Iter  t, const char * c, SizeType l){
    int base{16}; // The base, hexadecimal in this case
    // Use iterator_traits to deduce the type of the pointee.
    using value_type = typename std::iterator_traits<Iter>::value_type;
    value_type val;
    constexpr const std::size_t BYTE_SIZE{2ul}, VAL_SIZE{BYTE_SIZE*sizeof(decltype(val))};
    std::size_t length{static_cast<std::size_t>(l)};
    // Skip the hexadecimal prefix if any
    auto start=c;
    auto last=c + length;
    if constexpr ( DoCheck ) {
    	if( (c[0]=='0') && ((c[1]=='x') || (c[1]=='X')) ) start+=BYTE_SIZE;
    	// If length is not std::size_t, then is the size of the type to be converted, is not the buffer size.
    	if constexpr(!std::is_same<SizeType, std::size_t >::value ) {
    		last = start + length;
    	}
    }
    auto rval=std::tie(start,val);
    while( (start < last) && (start) ){
    	// Skip the test on hexadecimal prefix in hex_to_int, if already done
       	rval = hex_to_int<decltype(val),sizeof(decltype(val)), !DoCheck> ( start, val );
       	*t++=val;
    };

    return std::tuple<const char *, Iter>{start,t};
}

} // namespace string
} // namespace fcutils






#endif // HEX_TO_INT_HPP_ 
