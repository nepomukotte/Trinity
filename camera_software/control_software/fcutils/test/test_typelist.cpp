//
// test_typelist.cpp
//
//  Created on: 17 dic 2020
//      Author: cafagna
//


#include <iostream>
#include <type_traits>
#include <utility>
#include <tuple>
#include <variant>

#include "fcutils/type_traits/type_list.hpp"

namespace fct=fcutils::type_traits;

//template< typename Fun, typename Tuple, typename ... Arg, std::size_t ... Idxs>
//void apply_modifier_impl(Fun &&f, Tuple && t, std::index_sequence<Idxs ...>, Arg ... args ){
//        (...,f(std::get<Idxs>(t),args));
//}

template<typename T, T ... v > auto type_list_from_enum(){
//	enum: T {(v, ...)};
	return fct::type_list<std::integral_constant<T,v> ...>{};
//	return std::variant<std::integral_constant<T,v> ...>{};
}

template <typename Alternative, typename T> auto check(T t){
	if(t!=Alternative::value) {
		return 0ul;//" none ";
	} else {
		return static_cast<std::size_t>(Alternative::value);
	}
//	return (t!=Alternative::value)?0ul:static_cast<std::size_t>(Alternative::value);
}

template <typename Variant, typename T, std::size_t ... I > void test_loop_impl( Variant & v, T t, std::index_sequence<I ...>){
//	for( auto & v: {(std::get<I>(v), ... }){
//		(std::cout << ... << static_cast<std::size_t>(std::get<I>(v).value));
	(std::cout << ... << I );
		std::cout << std::endl;
		(std::cout << ... << static_cast<std::size_t>(std::variant_alternative_t<I,Variant>::value));
		std::cout << std::endl;
		(std::cout << ... << check<std::variant_alternative_t<I,Variant>>(t));
		std::cout << std::endl;
//	}
}

template <typename Variant, typename T, typename Indices = std::make_index_sequence<std::variant_size_v<Variant>>> void test_loop(Variant & v, T t){
	test_loop_impl(v,t,Indices{});
}
template <typename T> auto make_object( const T & v){

}
template<typename T, T ... v > auto create_list(){
//	enum: T {(v, ...)};
//	return fct::type_list<std::integral_constant<T,v> ...>{};
	return std::variant<std::integral_constant<T,v> ...>{};
}

template<typename T, T ... v > struct variant_from_enum{
	using type = std::variant<typename std::integral_constant<T,v>::type ...>;
};

template<typename T, T ... v > struct variant_integral_from_enum{
	using type = std::variant<std::integral_constant<T,v> ...>;
};

template<typename T, T ... v > auto create_type_list(){
	return fct::type_list<typename std::integral_constant<T,v>::type ...>{};
}

template <std::size_t Enum, typename TypeList> auto make_object_impl(std::index_sequence<Enum>, TypeList &){
	fct::get_type<static_cast<std::size_t>(Enum),TypeList> temp;
	return temp;
}

//template<typename T, T ... v > struct type_from_enum{
//	using typelist_type= fct::type_list<typename std::integral_constant<T,v>::type ...>;
//	using variant_type = std::variant<typename std::integral_constant<T,v>::type ...>;
//
//	auto enable_type()
//};
//



int main(){
	using tl=fct::type_list<int, double, float >;

	auto index=fct::list_search<bool,tl>::index;
	// Using the ntuple
	tl::type<2> f=4.5F;
	// Using the recursive metaprogramming
	fct::get_type<1,tl>::type d=3.14;

	// Uncomment next line to get an out of bound compile time error
//	fct::get_type<6,tl>::type junk;

	std::cout <<" Index:" << index <<"\n Type list struct size: " << sizeof(tl) << ", d: " << d << ", f: " << f
			<< " \n list_search size: " << sizeof(fct::list_search<double,tl>) << ", "<< sizeof(fct::list_search<bool,tl>)
			<< " \n get_type size: " << sizeof(fct::get_type<1,tl>) << ", "<< sizeof(fct::get_type<0,tl>)
			<< std::endl;

	enum class digit: std::size_t {one,two,three,four};
	auto tln=create_list<digit,digit::one,digit::two,digit::three,digit::four>();

//	auto tpl=std::make_variant(digit::one,digit::two,digit::three,digit::four);
	tln=std::integral_constant<digit,digit::three>{};
	std::visit([]( auto & v ){ std::cout << static_cast<std::size_t>(v()) << std::endl;}, tln);

	enum class my_type_flags: std::size_t {INT,FLOAT,DOUBLE};
	auto my_list=create_type_list<my_type_flags,my_type_flags::INT,my_type_flags::FLOAT,my_type_flags::DOUBLE>();
	variant_from_enum<my_type_flags,my_type_flags::INT,my_type_flags::FLOAT,my_type_flags::DOUBLE>::type my_variant;
	variant_integral_from_enum<my_type_flags,my_type_flags::INT,my_type_flags::FLOAT,my_type_flags::DOUBLE>::type my_integral_variant;

	test_loop(my_integral_variant,my_type_flags::FLOAT);

	return 0;
}




