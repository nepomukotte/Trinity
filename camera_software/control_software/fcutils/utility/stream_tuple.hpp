///
///
/// @brief To stream an tuple to the standard output
///
/// @file stream_tuple.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 3 feb 2022 
///

#ifndef UTILITY_STREAM_TUPLE_HPP_
#define UTILITY_STREAM_TUPLE_HPP_
#include <iosfwd> // forward declaration of the iostream objects
#include <tuple> // for std::tuple
#include <utility> // for std::index_sequence
#include <type_traits> // for is_pointer

// pretty-print a tuple

template<class Ch, class Tr, class Tuple, std::size_t... Is>
void stream_tuple_impl(std::basic_ostream<Ch,Tr>& os,
                      const Tuple& t,
                      std::index_sequence<Is...>)
{
	// Since C++17, folded expression
    ((os << (Is == 0? "" : ", ") << std::get<Is>(t)), ...);
}

template<class Ch, class Tr, class... Args>
auto& operator<<(std::basic_ostream<Ch, Tr>& os,
                 const std::tuple<Args...>& t)
{
    os << '(';
    stream_tuple_impl(os, t, std::index_sequence_for<Args...>{});
    return os << ')';
}




#endif // UTILITY_STREAM_TUPLE_HPP_ 
