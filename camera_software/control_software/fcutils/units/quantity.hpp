//!
//! quantity.h
//!
//!  Created on: 18 dic 2020
//!      Author: cafagna
//!
//! Types and utilities to set-up a basic and essential SI unit system.
//! It is based on the Bjarne Stroustrup The CPL, $28.7
//! Quantity is a value with an associated unit.
//!

#ifndef UNITS_QUANTITY_HPP_
#define UNITS_QUANTITY_HPP_

#include <units/units.hpp>
#include <string>
#include <iosfwd>

namespace fcutils{
namespace units{
namespace si{

template<typename U>
struct quantity {
	double val;
	explicit quantity(double d) : val{d} {}
	explicit operator double() const {return val;}
};

template<typename U>
const quantity<U> operator+(quantity<U> x, quantity<U> y){
	return quantity<U>{x.val+y.val};
}

template<typename U>
const quantity<U> operator-(quantity<U> x, quantity<U> y){
	return quantity<U>{x.val-y.val};
}

template<typename U>
const quantity<U> operator-(quantity<U> x){
	return quantity<U>{-x.val};
}

template<typename U>
const quantity<U> operator*(double x, quantity<U> y){
	return quantity<U>{x*y.val};
}

template<typename U>
const quantity<U> operator*(quantity<U> x, double y){
	return quantity<U>{y*x.val};
}

template<typename U>
const quantity<U> operator/(quantity<U> x, double y){
	return quantity<U>{x.val/y};
}

template<typename U1, typename U2>
const quantity<unit_plus<U1,U2>> operator*(quantity<U1> x, quantity<U2> y){
	return quantity<unit_plus<U1,U2>>{x.val*y.val};
}

template<typename U1, typename U2>
const quantity<unit_minus<U1,U2>> operator/(quantity<U1> x, quantity<U2> y){
	return quantity<unit_minus<U1,U2>>{x.val/y.val};
}

// User defined literals
// Scalar
quantity<scalar> operator"" _sc(long double d) { return quantity<scalar>{static_cast<double>(d)}; }
// meters
quantity<meter> operator"" _m(long double d) { return quantity<meter>{static_cast<double>(d)}; }
// seconds
quantity<second> operator"" _s(long double d) { return quantity<second>{static_cast<double>(d)}; }
// Ampere
quantity<ampere> operator"" _A(long double d) { return quantity<ampere>{static_cast<double>(d)}; }
// volts
quantity<volt> operator"" _V(long double d) { return quantity<volt>{static_cast<double>(d)}; }
// Kelvin
quantity<kelvin> operator"" _K(long double d) { return quantity<kelvin>{static_cast<double>(d)}; }
// Kelvin/second
quantity<unit_minus<kelvin,second>> operator"" _Kps(long double d) { return quantity<unit_minus<kelvin,second>>{static_cast<double>(d)}; }
// Volt/second
quantity<unit_minus<volt,second>> operator"" _Vps(long double d) { return quantity<unit_minus<volt,second>>{static_cast<double>(d)}; }
// Ampere/second
quantity<unit_minus<ampere,second>> operator"" _Aps(long double d) { return quantity<unit_minus<ampere,second>>{static_cast<double>(d)}; }

template<typename U> quantity<unit_plus<U,U>> square(quantity<U> x) { return quantity<unit_plus<U,U>>(x.val*x.val); }

template<typename U> bool operator==(quantity<U> x, quantity<U> y) { return x.val==y.val; }
template<typename U> bool operator!=(quantity<U> x, quantity<U> y) { return x.val!=y.val; }
template<typename U> bool operator<(quantity<U> x, quantity<U> y) { return x.val<y.val; }
template<typename U> bool operator<=(quantity<U> x, quantity<U> y) { return x.val<=y.val; }
template<typename U> bool operator>(quantity<U> x, quantity<U> y) { return x.val>y.val; }
template<typename U> bool operator>=(quantity<U> x, quantity<U> y) { return x.val>=y.val; }

template<typename U> constexpr bool is_negative(const quantity<U> &x) { return x.val < 0.; }

std::string suffix(int u, const char* x){
	std::string suf;
	if (u) {
		suf = ' ';
		suf += x;
		if (1<u) suf += '0'+u;
		if (u<0) {
			suf += '-';
			suf += '0'-u;
		}
	}
	return suf;
}

std::string suffix(int u, const char* x, bool lead_white_space ){
	std::string suf;
	if (u) {
		if(lead_white_space) suf = ' ';
		suf += x;
		if (1<u) suf += '0'+u;
		if (u<0) {
			suf += '-';
			suf += '0'-u;
		}
	}
	return suf;
}

template <typename U> std::string suffix(quantity<U> v){
	// rely on the overloaded suffix
	std::string suf{suffix(U::m,"m",false)};
	suf+=suffix(U::k,"K",false);
	suf+=suffix(U::a,"A",false);
	suf+=suffix(U::v,"V",false);
	suf+=suffix(U::s,"s",false);
	return suf;
}

template<typename U> std::ostream& operator<<(std::ostream& os, quantity<U> v) {
	return os << v.val << '[' <<suffix(U::m,"m") << suffix(U::k,"K") << suffix(U::a,"A")<< suffix(U::v,"V")<< suffix(U::s,"s")<<" ]";
}

std::ostream& operator<<(std::ostream& os, quantity<scalar> v) {
	return os << v.val ;
}

// Here we know how units are made.... do we need concepts?


} // namespace si
} // namespace units
} // namespace fcutils

#endif /* UNITS_QUANTITY_HPP_ */
