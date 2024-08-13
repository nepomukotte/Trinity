///
///
/// @brief A template class to just printout the statistics in the consumermanger. The Manager, must have a stats_print(std::ostream &) method.
/// NOTE: is not thread safe. pckts_ and processeds_ variables are not atomic in consumer_base
///
/// @file print_stats_consumer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 15 mar 2022 
///

#ifndef PRINT_STATS_CONSUMER_HPP_
#define PRINT_STATS_CONSUMER_HPP_

#include <functional>
#include <iosfwd>

#include "fcutils/consumer/consumer_base.hpp"

namespace fcutils{
namespace consumer{

// A template class to just printout the statistics in the consumermanger.
// The Manager, must have a stats_print(std::ostream &) method.
// NOTE: is not thread safe. pckts_ and processeds_ variables are not atomic in consumer_base
//
template <typename Data=unsigned char, template <typename D> class Checker = fcutils::checker::true_checker>
class print_stats_consumer: public consumer_base<Data,Checker>{
public:
	using self_type			=print_stats_consumer<Data,Checker>;
	using consumer_base_type=consumer_base<Data,Checker>;
	using checker_type 		=typename consumer_base_type::checker_type;
	using data_type			=typename consumer_base_type::data_type;
	using function_type		=std::function<void()>;

	template<typename Manager, typename C=checker_type>
	print_stats_consumer(const Manager & m, std::size_t frequency=100, std::ostream & os=std::cout, C && checker=checker_type{}):
		consumer_base_type{std::forward<checker_type>(checker)}, frequency_{frequency},
		function_{ [&manager=m, &o=os]() mutable { manager.print_stats(o); } } {};

	print_stats_consumer( print_stats_consumer && ) noexcept = default;
	self_type & operator=( print_stats_consumer && ) noexcept = default;

	void set_frequency(std::size_t frequency)  { frequency_=frequency;}

	self_type & operator()(const data_type * c, std::size_t length) override {
		consumer_base_type::update_stats(length);
		if( (consumer_base_type::pckts_%frequency_)==0 ) {
				function_();
		}
		return *this;
	}

private:
	function_type function_;
	std::size_t frequency_;
};

} // namespace consumer
} // namespace fcutils



#endif // PRINT_STATS_CONSUMER_HPP_ 
