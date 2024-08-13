///
///
/// @brief A new file
///
/// @file count_consumer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 15 mar 2022 
///

#ifndef COUNT_CONSUMER_HPP_
#define COUNT_CONSUMER_HPP_

#include <iosfwd>
#include <thread>
#include "fcutils/consumer/consumer_base.hpp"

namespace fcutils{
namespace consumer{

template <typename Data=unsigned char, template <typename D> class Checker = fcutils::checker::true_checker>
class count_consumer: public consumer_base<Data,Checker>{
public:
	using self_type			=count_consumer<Data,Checker>;
	using consumer_base_type=consumer_base<Data,Checker>;
	using checker_type 		=typename consumer_base_type::checker_type;
	using data_type			=typename consumer_base_type::data_type;

	template<typename C=checker_type>
	count_consumer(std::size_t frequency=100, C && checker=checker_type{}): consumer_base_type{std::forward<checker_type>(checker)}, frequency_(frequency) {};

	count_consumer( count_consumer && ) noexcept = default;
	self_type & operator=( count_consumer && ) noexcept = default;

	void set_frequency(std::size_t frequency)  { frequency_=frequency;}

	self_type &  operator()(const data_type * c, std::size_t length) override {
		consumer_base_type::update_stats(length);
		if( (consumer_base_type::pckts_%frequency_)==0 ) {
			std::cout << " Thread: " << std::this_thread::get_id() <<" . Packets processed: " << consumer_base_type::pckts_
					<< ", tot. bytes: " << consumer_base_type::processed_ << std::endl;
		}
		return *this;
	}

private:
	std::size_t frequency_;
};

} // namespace consumer
} //namespace fcutils



#endif // COUNT_CONSUMER_HPP_ 
