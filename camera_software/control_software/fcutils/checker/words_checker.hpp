///
///
/// @brief A new file
///
/// @file words_checker.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 29 ott 2021 
///

#ifndef WORDS_CHECKER_HPP_
#define WORDS_CHECKER_HPP_
#include "fcutils/checker/base_checker.hpp"
#include "fcutils/checker/return_policies.hpp"
#include "fcutils/type_traits/type_traits.hpp" // for is_bool
#include <array>
#include <utility> // for std::forward
#include <iostream>


namespace fcutils{
namespace checker{

///
/// A words checker. It does check for a word pattern.
/// The pattern must be present at the beginning of the word list.
/// the checker must be reset before using it again
/// The word pattern size is needed as a template parameter to leverage std::array usage.
/// The policy of the return value, is implemented via a return policy functor object.
/// The return policy functor object must implement an operator: Return operator()(bool,Data,std::size_t).
/// the checker operator(), does call the return policy using as actual parameter the result of its test, the data pointer (as passed by ASIO),
/// the counter of the data word tested.
///
template<std::size_t Size, typename Return=bool, template <typename, typename> class ReturnPolicy=do_nothing, typename Data=char, typename Key=char, bool Debug=false > struct words_checker final:
public base_checker<Data,Key,Return>{

	using self_type=words_checker<Size,Return,ReturnPolicy,Data,Key>;
	using base_type=base_checker<Data,Key,Return>;
	using array_type=std::array<typename base_type::data_type,Size>;
	using return_policy_type=ReturnPolicy<bool, typename base_type::return_type>;

	/// A constructor without the key. Just a list of parameter forwarded to the std::array constructor.
	template<typename ... Arg> words_checker(Arg&& ... args): test_pattern_{{std::forward<Arg>(args)...}} {
		static_assert( Size>=(sizeof...(Arg))," Ooopsss. List of parameters larger than std::array size !!!!!");
	};
	/// A constructor with the key and the return policy. Here we pass another array, instead of an initializer_list to enforce static check on array size
	/// We need to curly embrace the list to create the temporary std::array, i.e.: words_checker<3,char,char,bool> a{'k',{'a','b','c'}}
	words_checker(const typename base_type::key_type & k, const array_type & args, return_policy_type && rp=return_policy_type{}):
		base_type{k}, test_pattern_{args}, return_policy_{rp} {
	};

	// Copy
	words_checker(const self_type & c): base_type{c.key}, test_pattern_{c.test_pattern_}, return_policy_{c.return_policy_}{};
	self_type & operator=(const self_type & c){
		// Key is from the base type
		this->key=c.key;
		test_pattern_=c.test_pattern_;
		return_policy_=c.return_policy_;
	}

	typename base_type::return_type operator()( typename base_type::data_type * d) override {
		// Increment the counter of word read
		++word_count_;
		if constexpr (Debug) {
			std::cout << " Word checker. word_count_: " << word_count_ <<
					", passed: "<< passed_ <<
					", data: " << std::hex<< static_cast<unsigned int>(*d) << std::dec << std::endl;
		}
		// Check if we already found a good pattern
		if(passed_){
			// We do not need to check the test_pattern anymore, but we would like to check if the policy returns a zero, that is: read complete.
			// In case of read complete we need to reset the checker
			typename base_type::return_type to_be_read=return_policy_(passed_,d,word_count_);
			if(to_be_read==0) reset();
			return to_be_read;
		}
		// Here we haven't found a pattern but we have reached the end of the test_pattern
		// So we need to reset and start all over again ...
		if(iter_==test_pattern_.cend()) reset();
		// Assume test is not passed
		passed_=false;
		if(*d==*iter_){
			// a match found
			// it is the first word checked?
			if(iter_==test_pattern_.cbegin()) previous_=true;
			// if the previous check was fine, continue to consider previous_ true otherwise false
			previous_=(previous_)?previous_:false;
		}
		// Increment the test_pattern_ iterator
		++iter_;
		// If iterator reached the end then the full array have been tested
		if(iter_==test_pattern_.cend()) {
		// If previous_ is still true, then the test is passed otherwise is not passed
			passed_=previous_?true:false;
		}
		return return_policy_(passed_,d,word_count_);
	}

	typename base_type::return_type operator()( typename base_type::data_type * d, std::size_t l) override {
		typename base_type::return_type temp;
		for(std::size_t i=0; i<l;++i) {
			temp=operator()(d+i);
		}
		return temp;
	}

	void reset() override {
		passed_=previous_=false;
		iter_=test_pattern_.cbegin();
		word_count_=0;
	}

    std::size_t size(){ return test_pattern_.size(); }

private:
	const array_type test_pattern_;
	// Since C++17 begin and cbegin are constexpr
	typename array_type::const_iterator iter_{test_pattern_.cbegin()};
	std::size_t word_count_{};
	bool passed_{false},previous_{false};
	return_policy_type return_policy_{};
};

template<std::size_t Size> using bytes_checker=words_checker<Size,bool,do_nothing,char,char>;
template<std::size_t Size> using bytes_completer=words_checker<Size,std::size_t,do_nothing,char,char>;

} // namespace checker
} // namespace fcutils




#endif // WORDS_CHECKER_HPP_ 
