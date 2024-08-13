//
// test_chrono.cpp
//
//  Created on: 5 mag 2021
//      Author: cafagna
//

#include <chrono/timer.hpp>
#include <iostream>


namespace fcch=fcutils::chrono;

int main(){

	fcch::timer_nanosec tnano;
	fcch::timer_microsec tmicro;
	fcch::timer_millisec tmilli;
	fcch::timer_sec tsec;
	fcch::timer_minute tminute;
	fcch::timer_hour thour;
	fcch::timer_day tday;
	fcch::timer_week tweek;
	fcch::timer_month tmonth;
	fcch::timer_year tyear;

	std::string in;
	std::cout << " Enter anything ... (Ctrl-D to exit) \n";
	while(std::getline(std::cin,in)){
		tnano.print_diff(" Is a test, duration: ");
		tmicro.print_diff(" Is a test, duration: ");
		tmilli.print_diff(" Is a test, duration: ");
		tsec.print_diff(" Is a test, duration: ");
		tminute.print_diff(" Is a test, duration: ");
		thour.print_diff(" Is a test, duration: ");
		tday.print_diff(" Is a test, duration: ");
		std::cout << " A week duration: " << tweek;
		std::cout << " A month duration: " << tmonth;
		std::cout << " An year duration: " << tyear;
		std::cout << " Enter anything ... (Ctrl-D to exit) \n";
	}
	return 0;
}


