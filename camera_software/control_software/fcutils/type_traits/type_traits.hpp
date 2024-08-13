//
 // type_traits.h
 //
 //  Created on: 16 dic 2020
 //      Author: cafagna
 //

#ifndef TYPE_TRAITS_TYPE_TRAITS_HPP_
#define TYPE_TRAITS_TYPE_TRAITS_HPP_

#include <type_traits>
#include <iosfwd> // for std::(i)ostream forward declarations
#include <memory> // for std::shared_ptr

namespace fcutils{
namespace type_traits{

// Just to be sure that an output stream operator is defined for a class
// SFINAE in action here
template<typename T, typename S=std::ostream, typename = void>
struct is_to_stream_writable: std::false_type {};

template<typename T, typename S>
struct is_to_stream_writable<T, S,
        std::void_t< decltype( std::declval<S&>()<<std::declval<T>() )  > >
: std::true_type {};
// A touch of C++17 style ...
template<typename T, typename S=std::ostream> inline constexpr bool is_to_stream_writable_v = is_to_stream_writable<T,S>::value;

/// Just to be sure that an input stream operator is defined for a type
/// SFINAE in action here

template<typename T, typename S=std::istream, typename = void>
struct is_from_stream_readable: std::false_type {};

template<typename T, typename S>
struct is_from_stream_readable<T,S,
        std::void_t< decltype( std::declval<S&>()>>std::declval<T&>() )  > >
: std::true_type {};
// A touch of C++17 style ...
template<typename T, typename S> inline constexpr bool is_from_stream_readable_v = is_from_stream_readable<T,S>::value;

// Remove both CV and reference
template< typename T >
struct remove_cvref {
    using type=std::remove_cv_t< std::remove_reference_t<T> >;
};
// A touch of C++17 style ...
template < typename T> using remove_cvref_t = typename remove_cvref<T>::type;

// A check if a type is strictly boolean
template <typename T>
using is_bool = std::is_same< remove_cvref_t<T>, bool>;
template <typename T>
inline constexpr bool is_bool_v = is_bool<T>::value;

// Check if is a shared_ptr
// Forward declare shared_ptr
//template<typename T> class std::shared_ptr<T>;
template<typename T> struct is_shared_ptr : std::false_type {};
template<typename T> struct is_shared_ptr<std::shared_ptr<T>> : std::true_type {};
// A touch of C++17 style ...
template<typename T> inline constexpr bool is_shared_ptr_v = is_shared_ptr<T>::value;

// Is invocable with overload
template<typename T, typename Data, typename = void>
struct is_invocable_ovl: std::false_type {};

template<typename T, typename Data>
struct is_invocable_ovl<T, Data,
	std::void_t< decltype( std::declval<T>()(std::declval<const Data *>(),std::declval<std::size_t>()))>
	>
	: std::true_type {};
// A touch of C++17 style ...
template<typename T, typename Data> inline constexpr bool is_invocable_ovl_v = is_invocable_ovl<T,Data>::value;

// Has consume
template<typename T,  typename Data, typename = void>
struct has_consume: std::false_type {};

template<typename T, typename Data>
struct has_consume<T, Data,
	std::void_t< decltype( std::declval<T>().consume(std::declval<const Data *>(),std::declval<std::size_t>()))>
	>
	: std::true_type {};
// A touch of C++17 style ...
template<typename T, typename Data> inline constexpr bool has_consume_v = has_consume<T,Data>::value;

// Is explicit convertible.
// Many thanks to stackoverload for the code: https://stackoverflow.com/questions/16893992/check-if-type-can-be-explicitly-converted
template <class T, class U, class = void>
struct is_explicitly_convertible_to_impl : std::false_type {};

template <class T, class U>
struct is_explicitly_convertible_to_impl<
    T, U, std::void_t<decltype(static_cast<U>(std::declval<T>()))>>
    : std::true_type {};

template <class T, class U>
struct is_explicitly_convertible_to
    : is_explicitly_convertible_to_impl<T, U> {};

template <class T, class U>
inline constexpr bool is_explicitly_convertible_to_v =
    is_explicitly_convertible_to<T, U>::value;

/// const_compare. A struct for comparison to a constant.
// generic template, just declare a static with the value
template <typename T, T val, typename = void> struct const_compare{
        using value_type = T;
        using type = const_compare;
        constexpr static const value_type value = val;
        template <typename U>
        typename std::enable_if<is_explicitly_convertible_to_v<value_type,U>, bool>::type
        constexpr operator==(const U & u) const noexcept { return (static_cast<value_type>(u) == value); }
        constexpr operator value_type() const noexcept { return value; }
        constexpr value_type operator()() const noexcept { return value; }
};
// Specialization for enum case. In this case the enum underlying case is used as value_type.

template <typename T, T val>
struct const_compare<T,val,typename std::enable_if<std::is_enum<T>::value,std::void_t<T>>::type>{
        using value_type = T;
        using type = const_compare;
        using underlying_type = typename std::underlying_type<T>::type;
        enum: underlying_type {value = static_cast<underlying_type>(val)};
        template <typename U>
        typename std::enable_if<is_explicitly_convertible_to_v<underlying_type,U>, bool>::type
        constexpr operator==(const U & u) const noexcept { return (static_cast<underlying_type>(u) == value); }
        constexpr operator underlying_type() const noexcept { return static_cast<underlying_type>(value); }
        constexpr underlying_type operator()() const noexcept { return static_cast<underlying_type>(value); }
};

/// Cast an enumerator to the underlying type
template<typename E> // C++14
constexpr auto enum_to_underlying(E enumerator) noexcept
{
	return static_cast<std::underlying_type_t<E>>(enumerator);
}

} // namespace type_traits

} // namespace fcutils




#endif /* TYPE_TRAITS_TYPE_TRAITS_HPP_ */
