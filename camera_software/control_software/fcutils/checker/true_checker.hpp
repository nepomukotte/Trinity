///
///
/// @brief A new file
///
/// @file true_checker.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 28 ott 2021 
///

#ifndef TRUE_CHECKER_HPP_
#define TRUE_CHECKER_HPP_

#include "checker/base_checker.hpp"

namespace fcutils{
namespace checker{

// An always true checker. Does nothing.
template<typename Data=char, typename Key=char, typename Return=bool> struct true_checker final: public base_checker<Data,Key,Return>{
	using self_type=true_checker<Data,Key,Return>;
	using base_type=base_checker<Data,Key,Return>;
	// Let's just inherit all the base_checker creators
	using base_checker<Data,Key,Return>::base_checker;
	using base_checker<Data,Key,Return>::operator=;

	true_checker(true_checker && ) noexcept = default;
	self_type & operator=(true_checker &&) noexcept = default;


	typename base_type::return_type operator()( Data * ) override {
		// Implicit boolean conversion are used here ...
		return 1;
	}

	typename base_type::return_type operator()( Data *, std::size_t  ) override {
		// Implicit boolean conversion are used here ...
		return 1;
	}

	void reset() override {};
};


} // namespace checker
} // namespace fcutils



#endif // TRUE_CHECKER_HPP_ 
