///
///
/// @brief The base class for a consumer
///
/// @file consumer_base.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 25 feb 2022 
///

#ifndef CONSUMER_BASE_HPP_
#define CONSUMER_BASE_HPP_

#include <iostream>
#include <memory> // For std::shared_ptr
#include <thread> // For this_thread::get_id
#include <utility>

#include "fcutils/chrono/timer.hpp"
#include "fcutils/checker/base_checker.hpp"
#include "fcutils/checker/true_checker.hpp"


namespace fcutils{
namespace consumer{

template <typename Data=unsigned char, template<typename D> class Checker=fcutils::checker::true_checker> class consumer_base{
public:
	using duration_type	=std::chrono::milliseconds;
	using data_type		=Data;
	using data_ptr_type	=std::shared_ptr<data_type>;
	using checker_type	=Checker<data_type>;
	using self_type		=consumer_base;

	template <typename C=checker_type> consumer_base( C && checker=fcutils::checker::true_checker<data_type>{}): processed_{0}, pckts_{0},
			checker_{std::forward<checker_type>(checker)}{};
	virtual ~consumer_base(){};

	// Is a base class, so we need to get rid of the copy
	consumer_base(const consumer_base &) = delete;
	consumer_base & operator=(const consumer_base &) = delete;

	virtual self_type & operator()(const data_type * c, std::size_t  length) = 0;
	virtual void consume(const data_type * c, std::size_t length){
		this->operator()(c,length);
	}

	virtual void reset(){
		processed_=0;
		pckts_=0;
	}

	std::size_t processed(){ return processed_; }
	std::size_t packets(){ return pckts_; }

	virtual void print(std::ostream &o=std::cout) const {
		o<< " Consumer on thread :" << std::this_thread::get_id() << ", bytes processed: " << processed_ << ", number of packets: " << pckts_ << std::endl;
	}

	virtual void print_stats(std::ostream &os=std::cout, const std::string & prefix="") const {
		os<< prefix << " no. of processed packets: " << pckts_ << ", tot. Bytes: " << processed_ << '\n';
	}
protected:
	void update_stats(std::size_t length){
		processed_+=length;
		++pckts_;
	}
	// here some statistic
	std::size_t processed_;
	std::size_t pckts_;
	checker_type checker_;

};


} // namesapce consumer
} // namespace fcutils


#endif // CONSUMER_BASE_HPP_ 
