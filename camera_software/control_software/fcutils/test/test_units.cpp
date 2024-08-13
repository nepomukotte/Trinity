//
// test_units.cpp
//
//  Created on: 18 dic 2020
//      Author: cafagna
//
#include <units/quantity.hpp>
#include <iostream>
#include <chrono>


namespace fcusi=fcutils::units::si;
using fcusi::operator""_K;
using fcusi::operator""_s;


int main(){

	auto t=-12.3_K;
	auto d=1._s;
	std::chrono::duration<double> ds{3};
	fcusi::quantity<fcusi::second> dd{ds.count()};
	auto deriv=t/d;
	auto inc=deriv*static_cast<double>(std::chrono::steady_clock::now().time_since_epoch().count())*d;
    inc=deriv*dd;
    auto ratio=dd/d;
	std::cout << deriv << ", " << inc << ", " << ratio << std::endl;
	return 0;
}




