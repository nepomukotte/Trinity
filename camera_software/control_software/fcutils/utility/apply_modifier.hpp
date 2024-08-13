///
///
/// @brief An utility help function similar to the std::apply but invoking a function that modifies the std::tuple element. So, it returns void and leverage the folded expression
///
/// @file apply_modifier.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 3 feb 2022 
///

#ifndef UTILITY_APPLY_MODIFIER_HPP_
#define UTILITY_APPLY_MODIFIER_HPP_
#include <tuple> // for std::tuple_size, std::get
#include <utility> // for std::index, std::forward
#include <type_traits> // for std::remove_reference

namespace fcutils{

template< typename Fun, typename Tuple, typename ... Arg, std::size_t ... Idxs>
void apply_modifier_impl(Fun &&f, Tuple && t, std::index_sequence<Idxs ...>, Arg ... args ){
        (...,f(std::get<Idxs>(t),args));
}

template< typename Fun, typename Tuple, typename ... Arg>
void apply_modifier(Fun && f, Tuple && t, Arg && ... args){
        apply_modifier_impl(
        std::forward<Fun>(f), std::forward<Tuple>(t),
        std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>{},
        std::forward<Arg>(args) ... );
}


} // namespace fcutils


#endif // UTILITY_APPLY_MODIFIER_HPP_ 
