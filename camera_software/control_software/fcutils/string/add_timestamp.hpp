//
// add_timestamp.h
//
//  Created on: 19 dic 2020
//      Author: cafagna
//

#ifndef ADD_TIMESTAMP_HPP_
#define ADD_TIMESTAMP_HPP_
#include <string>
#include <string_view>
#include <chrono>
#include <ctime>

namespace fcutils {
namespace string {


enum class pos: bool { FRONT=true, BACK=false};
enum class conv: bool { CTIMEFMT=true, STRFTIME = false};

// convert time point with prefix to string:
template<fcutils::string::pos Position= fcutils::string::pos::FRONT, fcutils::string::conv Conversion= fcutils::string::conv::STRFTIME> std::string add_timestamp(std::string_view prefix,
		const std::chrono::system_clock::time_point &tp =
				std::chrono::system_clock::now()) {
	auto rawtime { std::chrono::system_clock::to_time_t(tp) };
	if constexpr (Conversion == fcutils::string::conv::CTIMEFMT){
		std::string temp;
		std::string_view ts;
		ts = std::ctime(&rawtime) ; // NOTE: not thread safe
		ts.remove_suffix(1); // skip trailing newline
		// Unfortunately no operator + available for string_view yet
		// Create a string and append the string_view to try to be a bit more efficient.
		temp.reserve(prefix.size() + ts.size());
		if constexpr(Position == fcutils::string::pos::FRONT) {
			temp = prefix;
			temp += ts;
		} else {
			temp = ts;
			temp += prefix;
		}
		return temp;
	} else {
		// A guess on the strftime size
		std::size_t strftimesize=20;
		std::string temp(prefix.size()+strftimesize,0);
		temp.resize(std::strftime(&temp[0], temp.size(),
				"%Y_%b_%d_%H_%M_%S", std::localtime(&rawtime)));
		if constexpr(Position == fcutils::string::pos::FRONT) {
			temp.insert(0,prefix);
		} else {
			temp += prefix;
		}
		return temp;
	}
}

} // namespace string
} // namespace fcutils
#endif /* ADD_TIMESTAMP_HPP_ */
