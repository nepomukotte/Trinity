///
///
/// @brief A new file
///
/// @file progam_name_and_timestamp.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 21 mar 2022 
///

#ifndef PROGRAM_NAME_AND_DATE_HPP_
#define PROGRAM_NAME_AND_DATE_HPP_

#include <string>
#include <chrono>

#include "fcutils/string/add_timestamp.hpp"

namespace fcutils{
namespace string{

class program_name_and_date{
public:
	explicit program_name_and_date(const std::string & c=program_invocation_short_name, const std::string & post=" -->> "): progname_{c}, post_{post} {
		sanity_check();
	};

	program_name_and_date(const program_name_and_date & ) = default;
	program_name_and_date& operator=(const program_name_and_date & ) = default;
	program_name_and_date(program_name_and_date && ) noexcept = default;
	program_name_and_date& operator=(program_name_and_date && ) noexcept = default;

	const char * c_str() const {
		stemp_ = fcutils::string::add_timestamp<fcutils::string::pos::FRONT, fcutils::string::conv::STRFTIME> (progname_) + post_;
		return stemp_.c_str();
	}

	void set_name(const std::string & c){ progname_=c; sanity_check();}
	void set_postfix(const std::string & c){ post_=c;}

	const char * get_name() const { return progname_.c_str(); }

	operator std::string() const {
		return std::string{c_str()};
	}

private:
	std::string progname_;
	mutable std::string stemp_;
	std::string post_;

	void sanity_check(){
		// Check if the string progname_ starts and ends with a whitespace. If no, add it.
		if(!progname_.empty()){
			const char WS=' ';
			if( progname_[0]!=WS ){
				progname_.insert(0,1,WS);
			}
			if( progname_[progname_.size()-1]!=WS ){
				progname_+=WS;
			}
		}
	}
};

std::ostream & operator << (std::ostream & o, program_name_and_date & p){
	o<< p.c_str();
	return o;
}

} // namespace string
} // namespace fcutils




#endif // PROGRAM_NAME_AND_DATE_HPP_ 
