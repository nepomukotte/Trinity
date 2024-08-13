//
// flat_random.h
//
//  Created on: 16 dic 2020
//      Author: cafagna
//

#ifndef RANDOM_FLAT_RANDOM_HPP_
#define RANDOM_FLAT_RANDOM_HPP_

#include <random> // random_engine and distribution
#include <chrono> // for clock facilities
#include <type_traits> // for std::make_unsigned
#include <cstdint> // for std::uint32_t

namespace fcutils{
namespace random{

template <typename ReturnType=int, template <typename> class Distribution=std::uniform_int_distribution, typename Engine=std::default_random_engine > class rand_flat{
public:
	using return_type=ReturnType;
	using distribution_type=Distribution<ReturnType>;
	using engine_type=Engine;
	using seed_type=std::uint32_t;
	using parameter_type=typename distribution_type::param_type;
	// Creators
	rand_flat(const return_type & low, const return_type & high, const seed_type & seed =
			static_cast<
				std::make_unsigned<decltype(std::chrono::steady_clock::now().time_since_epoch().count())>::type
			> (std::chrono::steady_clock::now().time_since_epoch().count())
	): dist_{low,high}, re_{seed} { };

	rand_flat(const parameter_type & par, const seed_type & seed =
			static_cast<
				std::make_unsigned<decltype(std::chrono::steady_clock::now().time_since_epoch().count())>::type
			> (std::chrono::steady_clock::now().time_since_epoch().count())
	): dist_{par}, re_{seed} { };

	// Creators
	template <typename SSequence> rand_flat(const return_type & low, const return_type & high, SSequence & seeds): dist_{low,high}, re_{seeds} { };

	template <typename SSequence> rand_flat(const parameter_type & par, SSequence & seeds): dist_{par}, re_{seeds} { };
	// Use the operator () to return a random
	// number in the range [low,high]
	return_type operator()(){ return dist_(re_);}
	void seed(seed_type s){ re_.seed(s);}
	template <typename SSequence> void seed(SSequence &s){ re_.seed(s); }
	void discard( unsigned long long ncycle ){ re_.discard(ncycle); }
	return_type get_low() const { return dist_.a(); }
	return_type get_high() const { return dist_.b(); }
	parameter_type get_params() const { return dist_.param(); }
private:
	engine_type re_;
	distribution_type dist_;
}; // class rand_flat

using rand_int=rand_flat<int,std::uniform_int_distribution,std::default_random_engine>;
using rand_double=rand_flat<double,std::uniform_real_distribution,std::default_random_engine>;


} // namespace random
} // namespace fcutils


#endif /* RANDOM_FLAT_RANDOM_HPP_ */
