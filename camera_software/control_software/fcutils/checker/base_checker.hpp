///
///
/// @brief A new file
///
/// @file base_checker.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 28 ott 2021 
///

#ifndef BASE_CHECKER_HPP_
#define BASE_CHECKER_HPP_
#include <utility> // for std::move
#include <type_traits> // for std::remove_cv and reference

namespace fcutils{
namespace checker{

// A base checker for a collection of checkers.
// Data is the type of data to be checked,
// Key is the type of key to be used in case an associative containers is needed
template<class Data=char, class Key=char, typename Return=bool> struct base_checker{
	// Remove both CV and reference C++11 style ...
	// We would like to be sure that Data an Key are just the "bold" type.
	using data_type=typename std::remove_cv< typename std::remove_reference<Data>::type >::type;
	using key_type=typename std::remove_cv< typename std::remove_reference<Key>::type >::type;
	using return_type=Return;
	// As usual remove the copy ... is a base class
	base_checker(const base_checker &) = delete;
	base_checker & operator=(const base_checker &) = delete;
	// ... and make the destructor virtual
	virtual ~base_checker(){};
	// Because we are declaring ctor or dtor ... best create the standard constructor ...
	base_checker()= default;
	// At this point ... let's add a constructor with the key
	base_checker(key_type && k ): key{std::move(k)}{};
	base_checker(const key_type & k): key{k}{};
	// ... and move constructors
	base_checker(base_checker && bc) noexcept : key{std::move(bc.key)}{};
	base_checker & operator=(base_checker && bc) noexcept {
		key=std::move(bc.key);
		return *this;
	}
	// Type erasure is here ... a classical () operator will do the job.
	virtual return_type operator()(data_type *)=0;
	virtual return_type operator()(data_type *,std::size_t)=0;
	virtual void reset()=0;
	// Setter
	void set_key( const key_type & k ){
		key=k;
	}
	// Data member ... eventually
	key_type key{};
};

} // namespace checker
} // namespace fcutils

#endif // BASE_CHECKER_HPP_ 
