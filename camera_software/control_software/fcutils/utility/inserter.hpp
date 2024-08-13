///
///
/// @brief To "insert" a value inside an integer given an offset and a mask
///
/// @file inserter.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 15 feb 2022 
///

#ifndef UTILITY_INSERTER_HPP_
#define UTILITY_INSERTER_HPP_
#include <utility>

namespace fcutils{

template<typename Target, typename Source, typename Mask> constexpr Target insert(Target & t, Source s, int offset, Mask m){
	s=s&static_cast<Source>(m);
	s=s<<offset;
	return t=static_cast<Target>(s)|t;
}

template<typename Target, typename Source, typename Mask> constexpr Target insert(Target * t, Source s, int offset, Mask m){
	s=s&static_cast<Source>(m);
	s=s<<offset;
	return *t=static_cast<Target>(s)|*t;
}

// Forward declare the fcutils mask type, just to avoid dependence from integer.hpp
// The mask type must declare an unit_type as the underlying type of the actual mask
template <int LastBit,int Offset, int BitSize>
struct mask_t;

template<int BitSize, int MaskSize, int Offset, template<int, int, int> class Mask=fcutils::mask_t>  struct inserter{
	// Mask type must declare the mask type as uint_type
	using mask_type=Mask<MaskSize,0,BitSize>;
	using uint_type=typename mask_type::uint_type;
	enum: uint_type {mask=static_cast<uint_type>(mask_type::mask)};

	template <typename Target, typename Source> Target insert(Target & t, Source && s) const {
		return fcutils::insert<Target,Source,uint_type>(t,std::forward<Source>(s),Offset,mask);
	}

	template <typename Target, typename Source> Target insert(Target * t, Source && s) const {
		return fcutils::insert<Target,Source,uint_type>(t,std::forward<Source>(s),Offset,mask);

	}

	uint_type get_mask() const { return static_cast<uint_type>(mask); }
};

} // namespace fcutils




#endif // UTILITY_INSERTER_HPP_ 
