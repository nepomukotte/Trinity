//
 // clamp.hpp
 //
 //  Created on: 18 dic 2020
 //      Author: cafagna
 //

#ifndef CLAMP_HPP_
#define CLAMP_HPP_
/// The idea and interface is based on the Boost.Algorithm library
/// authored by Marshall Clow.
/// This is just implemented to avoid the usage of C++17 std:clamp. No standard support for this ISO release yet.

#include <algorithm> // for std::less
#include <cassert> // for run-time assert macro

namespace fcutils{
namespace algorithm {

template<class T>
constexpr const auto & clamp( const T& val, const T& lo, const T& hi) {
	assert( !(hi < lo) );
	return (val < lo) ? lo : (hi < val) ? hi : val;
}

template<class T, class Compare=std::less<T> >
constexpr const auto & clamp( const T& val, const T& lo, const T& hi, Compare comp )
{
	assert( !comp(hi, lo) );
	return comp(val, lo) ? lo : ( comp(hi, val) ? hi : val) ;
}

} // namespace algorithm
}  // namespace fcutils





#endif /* CLAMP_HPP_ */
