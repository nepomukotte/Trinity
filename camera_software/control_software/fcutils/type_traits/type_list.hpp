//
 // type_list.h
 //
 //  Created on: 16 dic 2020
 //      Author: cafagna
 //

#ifndef TYPE_TRAITS_TYPE_LIST_HPP_
#define TYPE_TRAITS_TYPE_LIST_HPP_
#include <tuple> // for std::tuple
#include <cstddef> // for std::size_t

namespace fcutils{
namespace type_traits{

// A simple type list to automatic number the type
// It is extracted from Alexandrescu's book and updated to C++17
// Cppannotation helped: https://gitlab.com/fbb-git/cppannotations
//
template <class... Args>
struct type_list
{
	constexpr type_list() = default;
	type_list(type_list const &) = delete;
	enum { size = sizeof ...(Args) };
	template <std::size_t N>
	using type = typename std::tuple_element<N, std::tuple<Args...>>::type;
};

// An alternative to the ntuple usage would use recursivity
// Primary
template<std::size_t index, typename TypeList>
struct get_type;

// Specialization
// Type not found or empty list
template <std::size_t index>
struct get_type<index, type_list<>>
{
	static_assert(index < 0, " get_type index out of bounds.");
	using type=get_type;
};

// Index == 0 case
template <typename Head, typename ...Tail>
struct get_type<0, type_list<Head, Tail...>>
{
	using type=Head;
};

// Other index values
template <std::size_t index, typename Head, typename ...Tail>
struct get_type<index, type_list<Head, Tail...>>
{
	using type= typename get_type<index - 1, type_list<Tail...>>::type;
};

// Primary
template <typename ...Types>
struct list_search
{
	list_search(list_search const &) = delete;
};

// Specialization
// Empty list
template <typename SearchType>
struct list_search<SearchType, type_list<> >
{
	list_search(list_search const &) = delete;
	enum { index = -1 };
};

// First type
template <typename SearchType, typename ... Tail>
struct list_search< SearchType, type_list<SearchType, Tail...> >
{
	list_search(list_search const &) = delete;
	enum { index = 0};
};
// Get the index

template <typename SearchType, typename Head, typename ...Tail>
struct list_search<SearchType, type_list<Head, Tail...> >
{
	list_search(list_search const &) = delete;
private:
	enum {tmp = list_search<SearchType, type_list<Tail...>>::index};
public:
	enum {index = tmp == -1 ? -1 : 1 + tmp};
};

// We need to add some more facilities to manage type list ... like purge duplicate, add or remove types etc. etc
} // namespace type_traits

} // namespace fcutils


#endif /* TYPE_TRAITS_TYPE_LIST_HPP_ */
