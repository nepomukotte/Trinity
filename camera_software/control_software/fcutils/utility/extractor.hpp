///
///
/// @brief To "extract" an integer with a given bit size and bit offset
///
/// @file extractor.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 14 feb 2022 
///

#ifndef UTILITY_EXTRACTOR_HPP_
#define UTILITY_EXTRACTOR_HPP_


namespace fcutils{

template<typename Source, typename Mask> constexpr Source extract(const Source & s, int offset, const Mask & m){
	Source t{s};
	t=t>>offset;
	return t&m;
}

template<typename Source, typename Mask> constexpr Source extract(const Source * s, int offset, const Mask & m){
	Source t{*s};
	t=t>>offset;
	return t&m;
}

// Forward declare the fcutils mask type, just to avoid dependence from integer.hpp
// The mask type must declare an unit_type as the underlying type of the actual mask
template <int LastBit,int Offset, int BitSize>
struct mask_t;

template<int BitSize, int MaskSize, int Offset, template<int, int, int> class Mask=fcutils::mask_t>  struct extractor{
	// Mask type must declare the mask type as uint_type
	using mask_type=Mask<MaskSize,0,BitSize>;
	using uint_type=typename mask_type::uint_type;
	enum: uint_type {mask=static_cast<uint_type>(mask_type::mask)};

	template <typename Source> Source extract(const Source & s) const {
		return fcutils::extract<Source,uint_type>(s,Offset,mask);
	}

	template <typename Source> Source extract(const Source * s) const {
		return fcutils::extract<Source,uint_type>(s,Offset,mask);
	}

	uint_type get_mask() const { return static_cast<uint_type>(mask); }
};

} // namespace fcutils

#endif // UTILITY_EXTRACTOR_HPP_ 
