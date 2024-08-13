//
// timer.h
//
//  Created on: 5 mag 2021
//      Author: cafagna
//

#ifndef CHRONO_TIMER_HPP_
#define CHRONO_TIMER_HPP_

// Inspired by N. Jousuits

#include <chrono> // for time points and durations
#include <string>
#include <iostream>
#include <type_traits> // for std::is_void_v

// Just an unnamed namespace to limit the alias to the current unit
namespace {
namespace ch = std::chrono;
// namespace alias
}

namespace fcutils {
namespace chrono {

template <typename Duration> std::string suffix(){
	std::string suf{"?????"};
	if constexpr (std::is_same<typename Duration::period,std::nano>::value) suf="nsec";
	if constexpr (std::is_same<typename Duration::period,std::micro>::value) suf="usec";
	if constexpr (std::is_same<typename Duration::period,std::milli>::value) suf="msec";
	if constexpr (std::is_same<typename Duration::period,std::ratio<1>>::value) suf="sec";
	if constexpr (std::is_same<typename Duration::period,std::ratio<60>>::value) suf="min";
	if constexpr (std::is_same<typename Duration::period,std::ratio<3600>>::value) suf="hour";
	if constexpr (std::is_same<typename Duration::period,std::ratio<86400>>::value) suf="day";
	if constexpr (std::is_same<typename Duration::period,std::ratio<604800>>::value) suf="week";
	if constexpr (std::is_same<typename Duration::period,std::ratio<2629746>>::value) suf="month";
	if constexpr (std::is_same<typename Duration::period,std::ratio<31556952>>::value) suf="year";
	return suf;
}

// Unfortunately no string litteral can be passed as a template value... until C++20, so we will use an external function.
template<typename Time_difference, typename Clock=ch::high_resolution_clock > class timer {
public:
	using time_diff_type=Time_difference;
	using clock_type=Clock;
	using time_point_type= typename clock_type::time_point;

	timer() :
			last_ { clock_type::now() } {
	}
	;
	void reset() {
		last_ = clock_type::now();
	}
	auto diff() const {
		auto now { clock_type::now() };
		auto time_diff { (now - last_) };
		return ch::duration_cast<time_diff_type>(time_diff);
	}

	auto diff_reset() {
		auto time_diff=diff();
		reset();
		return ch::duration_cast<time_diff_type>(time_diff);
	}

	void print_diff( const std::string &msg = "Time difference: ", std::ostream &o = std::cout) const {
		auto time_diff = diff();
		o << msg << time_diff.count() << suffix<time_diff_type>() << '\n';
		return;
	}

private:
	time_point_type last_;
};


using timer_nanosec=timer<ch::nanoseconds>;
using timer_microsec=timer<ch::microseconds>;
using timer_millisec=timer<ch::milliseconds>;
using timer_sec=timer<ch::seconds>;
using timer_minute=timer<ch::minutes>;
using timer_hour=timer<ch::hours>;
// Anticipating a bit C++20 ...
// Using a double to get a meaningful value for short durations.
using timer_day=timer<ch::duration<double,std::ratio<86400>>>;
using timer_week=timer<ch::duration<double,std::ratio<604800>>>;
using timer_month=timer<ch::duration<double,std::ratio<2629746>>>;
using timer_year=timer<ch::duration<double,std::ratio<31556952>>>;

template <class CharT, class Traits, class Duration>

std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits> & o, const timer<Duration> & t ){
			auto time_diff = t.diff();
			return o << time_diff.count() << suffix<Duration>() << '\n';
}

} // namespace chrono
} // namesapce fcutils

#endif /* CHRONO_TIMER_HPP_ */
