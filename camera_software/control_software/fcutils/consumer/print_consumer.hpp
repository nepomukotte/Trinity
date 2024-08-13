///
///
/// @brief A new file
///
/// @file print_consumer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 15 mar 2022 
///

#ifndef PRINT_CONSUMER_HPP_
#define PRINT_CONSUMER_HPP_

#include <iostream>
#include <iterator>
#include <algorithm>
#include <thread>
#include <utility>

#include "fcutils/consumer/consumer_base.hpp"

namespace fcutils{
namespace consumer{

enum class PRINT_FORMAT: bool { BINARY=true, TEXT=false };

template <typename Data=unsigned char, template <typename D> class Checker = fcutils::checker::true_checker, PRINT_FORMAT Format=PRINT_FORMAT::TEXT>
struct print_consumer: public consumer_base<Data,Checker>{
	using self_type			=print_consumer<Data,Checker,Format>;
	using consumer_base_type=consumer_base<Data,Checker>;
	using checker_type 		=typename consumer_base_type::checker_type;
	using data_type			=typename consumer_base_type::data_type;

	using consumer_base_type::consumer_base;
	using consumer_base_type::operator=;

	print_consumer( print_consumer && ) noexcept = default;
	self_type & operator=( print_consumer && ) noexcept = default;

	self_type & operator() (const data_type * c, std::size_t length) override {
		std::cout << " Thread: " << std::this_thread::get_id() <<". Consuming a message of size (bytes): " << length <<", values:\n";
		if constexpr ( Format == PRINT_FORMAT::TEXT ){
			std::ostreambuf_iterator<data_type> oo (std::cout.rdbuf()); // output iterator for cout
			std::copy(c,c+length,oo);
		} else {
			for(int i=0; i<length; ++i) std::cout << std::hex<< std::to_integer<short>(static_cast<std::byte>(static_cast<unsigned char>(*(c+i)))) <<' ';
		}
		std::cout << std::dec << std::endl;
		consumer_base_type::update_stats(length);
		return *this;
	}
};

template<typename Data=unsigned char, template <typename D> class Checker=fcutils::checker::true_checker>
using printbinary_consumer = print_consumer<Data,Checker,PRINT_FORMAT::BINARY>;
template<typename Data=unsigned char, template <typename D> class Checker=fcutils::checker::true_checker>
using printtext_consumer = print_consumer<Data,Checker,PRINT_FORMAT::TEXT>;

} // namespace fcutils
} // namesapce consumer




#endif // PRINT_CONSUMER_HPP_ 
