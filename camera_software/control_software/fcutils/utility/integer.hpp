///
///
/// @brief An helper class to get the smaller integer with at least the number of bits in the template
///
/// @file integer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 14 feb 2022 
///

#ifndef UTILITY_INTEGER_HPP_
#define UTILITY_INTEGER_HPP_

#include <cstdint> // for intX_t, uintX_t, and fast and least version

namespace fcutils{

// General template
template <int Bits>
  struct signed_base_helper{};
//Specializations
template <> struct signed_base_helper<8> {
	typedef int8_t exact;
	typedef int_fast8_t fast;
	typedef int_least8_t least;
};
template <> struct signed_base_helper<16> {
	typedef int16_t exact;
	typedef int_fast16_t fast;
	typedef int_least16_t least;
};
template <> struct signed_base_helper<32> {
	typedef int32_t exact;
	typedef int_fast32_t fast;
	typedef int_least32_t least;
};
template <> struct signed_base_helper<64> {
	typedef int64_t exact;
	typedef int_fast64_t fast;
	typedef int_least64_t least;
};

// General template
template <int Bits>
  struct unsigned_base_helper{};
//Specializations
template <> struct unsigned_base_helper<8> {
	typedef uint8_t exact;
	typedef uint_fast8_t fast;
	typedef uint_least8_t least;
};
template <> struct unsigned_base_helper<16> {
	typedef uint16_t exact;
	typedef uint_fast16_t fast;
	typedef uint_least16_t least;
};
template <> struct unsigned_base_helper<32> {
	typedef uint32_t exact;
	typedef uint_fast32_t fast;
	typedef uint_least32_t least;
};
template <> struct unsigned_base_helper<64> {
	typedef uint64_t exact;
	typedef uint_fast64_t fast;
	typedef uint_least64_t least;
};

// Using C++17
template <int Bits> constexpr int min_bit(){
	enum: int {EIGHT=8,SIXTEEN=16,THIRTYTWO=32,SIXTYFOUR=64};
	static_assert((Bits<=SIXTYFOUR),"Error in fcutils::min_bit, number of bits larger than 64 !!!!");
	if constexpr (Bits<=EIGHT){
		return static_cast<int>(EIGHT);
	} else {
		if constexpr (Bits <= SIXTEEN ) {
			return static_cast<int>(SIXTEEN);
		} else {
			if constexpr (Bits <= THIRTYTWO ){
				return static_cast<int>(THIRTYTWO);
			} else {
				return static_cast<int>(SIXTYFOUR);
			}
		}
	}
}

// Using C++11
constexpr int min_bit(const int & b){
	enum: int {EIGHT=8,SIXTEEN=16,THIRTYTWO=32,SIXTYFOUR=64};
	auto val=SIXTYFOUR;
	if(b<=THIRTYTWO) val=THIRTYTWO;
	if(b<=SIXTEEN) val=SIXTEEN;
	if(b<=EIGHT) val=EIGHT;
	return static_cast<int>(val);
}

// Using C++17
//template<int Bits>
//struct int_t: signed_base_helper<min_bit<Bits>()> {};
//
//template<int Bits>
//struct uint_t: unsigned_base_helper<min_bit<Bits>()> {};

// Using C++11
template<int Bits>
struct int_t: signed_base_helper<min_bit(Bits)> {};

template<int Bits>
struct uint_t: unsigned_base_helper<min_bit(Bits)> {};

template<int LastBit, int Offset=0, int BitSize=LastBit, typename I=typename uint_t<BitSize>::least> constexpr I mask(int b=LastBit-1) {
	return (b>Offset)?(I{1u}<<b)+mask<LastBit,Offset,BitSize,I>(b-1):(I{1u}<<Offset);
}

template<int BitSize, template<int> class I=fcutils::uint_t> constexpr typename I<BitSize>::least mask(int lastbit=BitSize-1,int offset=0) {
	using uint_type= typename I<BitSize>::least;
	return (lastbit>offset)?(uint_type{1u}<<lastbit)+mask<BitSize,I>(lastbit-1,offset):(uint_type{1u}<<offset);
}

template <int LastBit,int Offset=0, int BitSize=LastBit>
struct mask_t {
	using uint_type=typename uint_t<BitSize>::least;
	enum: uint_type {mask=mask<LastBit,Offset,BitSize,uint_type>()};
	operator uint_type() const { return static_cast<uint_type>(mask); }
	// Force to explicit to avoid implicit conversion to any type ....
	template<typename Cast> explicit operator Cast() const { return static_cast<Cast>(static_cast<uint_type>(mask)); }
};


} // namespace fcutils




#endif // UTILITY_INTEGER_HPP_ 
