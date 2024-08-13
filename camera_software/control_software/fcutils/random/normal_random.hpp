//
// normal_random.h
//
//  Created on: 16 dic 2020
//      Author: cafagna
//

#ifndef RANDOM_NORMAL_RANDOM_HPP_
#define RANDOM_NORMAL_RANDOM_HPP_

#include <random> // for std::default_random_engine and std::normal_distribution
#include <chrono> // for clock facilities
#include <type_traits> // form std::make_unsigned
#include <cstdint> // for std::uint32_t

namespace fcutils{
namespace random{

template <typename ReturnType=double, template <typename> class Distribution=std::normal_distribution, typename Engine=std::default_random_engine > class grand{
public:
	using return_type=ReturnType;
	using distribution_type=Distribution<ReturnType>;
	using engine_type=Engine;
	using seed_type=std::uint32_t;
	using parameter_type=typename distribution_type::param_type;
	// Creators
	grand(const return_type & mean, const return_type & sigma, const seed_type & seed =
			static_cast<
				std::make_unsigned<decltype(std::chrono::steady_clock::now().time_since_epoch().count())>::type
			> (std::chrono::steady_clock::now().time_since_epoch().count())
	): dist_{mean,sigma}, re_{seed} { };

	grand(const parameter_type & par, const seed_type & seed =
			static_cast<
				std::make_unsigned<decltype(std::chrono::steady_clock::now().time_since_epoch().count())>::type
			> (std::chrono::steady_clock::now().time_since_epoch().count())
	): dist_{par}, re_{seed} { };

	// Creators
	template <typename SSequence> grand(const return_type & mean, const return_type & sigma, SSequence & seeds): dist_{mean,sigma}, re_{seeds} { };

	template <typename SSequence> grand(const parameter_type & par, SSequence & seeds): dist_{par}, re_{seeds} { };
	// Use the operator () to return a random
	// number in the range [low,high]
	return_type operator()(){ return dist_(re_);}
	void seed(seed_type s){ re_.seed(s);}
	template <typename SSequence> void seed(SSequence &s){ re_.seed(s); }
	void discard( unsigned long long ncycle ){ re_.discard(ncycle); }
	return_type get_mean() const { return dist_.mean(); }
	return_type get_sigma() const { return dist_.stddev(); }
	parameter_type get_params() const { return dist_.param(); }
private:
	std::default_random_engine re_;
	distribution_type dist_;
}; // class grand




} // namespace random
} // namespace fcutils


#endif /* RANDOM_NORMAL_RANDOM_HPP_ */
