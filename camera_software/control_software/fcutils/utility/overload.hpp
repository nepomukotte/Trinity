///
///
/// @brief The overload pattern. Very usefull in case you have to visit a std::variant
///
/// @file oveload.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 16 dic 2021 
///

#ifndef UTILITY_OVERLOAD_HPP_
#define UTILITY_OVERLOAD_HPP_

namespace fcutils{

// ”inherit” all function call operators of passed base types:
template<typename... Base>
struct overload : Base...
{
	using Base::operator()...;
};

// Deduction guides. Since C++17. Needed to auto deduce types for a class.
// base types are deduced from passed arguments:
template<typename... Base>
overload(Base...) -> overload<Base...>;

} // namespace fcutils

#endif // UTILITY_OVERLOAD_HPP_ 
