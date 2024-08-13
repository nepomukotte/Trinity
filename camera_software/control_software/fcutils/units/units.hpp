//!
//! units.h
//!
//!  Created on: 18 dic 2020
//!      Author: cafagna
//!
//! Types and utilities to set-up a basic and essential SI unit system.
//! It is based on the Bjarne Stroustrup The CPL, $28.7
//! Units can be devided, then minus sign is used, or added, then the plus is used.
//!

#ifndef UNITS_UNITS_HPP_
#define UNITS_UNITS_HPP_

#include <iostream>

namespace fcutils {
namespace units {
namespace si {

template<int Meter, int Second, int Ampere, int Volt, int Kelvin>
struct unit {
	enum {
		m = Meter, k = Kelvin, a = Ampere, v = Volt,  s = Second
	};
};

// A couple of type functions to get the right dimension in case of operation between units.
// Units should be summed up in case of multiplication ...
template<typename U1, typename U2>
struct uplus {
	using type= unit<U1::m+U2::m ,U1::s+U2::s, U1::a+U2::a, U1::v+U2::v, U1::k+U2::k>;
};
template<typename U1, typename U2> using unit_plus = typename uplus<U1,U2>::type;

// ... and subtracted in case of division.
template<typename U1, typename U2>
struct uminus {
	using type= unit<U1::m-U2::m ,U1::s-U2::s, U1::a-U2::a, U1::v-U2::v, U1::k-U2::k>;
};
template<typename U1, typename U2> using unit_minus = typename uminus<U1,U2>::type;


using scalar=unit<0,0,0,0,0>;
// Basic units
using meter= unit<1,0,0,0,0>;
using second=unit<0,1,0,0,0>;
using ampere=unit<0,0,1,0,0>;
using volt=  unit<0,0,0,1,0>;
using kelvin=unit<0,0,0,0,1>;
// Composed units
// meter/second
using mps =unit<1,-1,0,0,0>;
// ampere/second
using aps =unit<0,-1,1,0,0>;
// voltage/second
using vps =unit<0,-1,0,1,0>;
// kelvin/second
using kps =unit<0,-1,0,0,1>;

} // namespace si

} // namespace units

} // namespace fcutils



#endif /* UNITS_UNITS_HPP_ */
