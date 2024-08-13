///
///
/// @brief A base class for header or trailers. Constructed as a derived class of a std::array
///
/// @file header_trailer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 8 mar 2022 
///

#ifndef BASE_HT_HPP_
#define BASE_HT_HPP_
#include <array>
#include <utility> // For std::forward

namespace fcutils{
namespace message{

template <typename Data, std::size_t Size> struct base_ht{
	using data_type		=Data;
	using data_ptr_type	=Data *;
	using self_type		=base_ht<Data,Size>;
	using array_type	=std::array<Data,Size>;
	// std::array constructors and assignments are not defined because it follows the rule of an aggregate type. So we would use a little trick ....
	template <typename ... Arg,
			typename std::enable_if<sizeof...(Arg) <= Size>::type * = nullptr>
		   constexpr base_ht (Arg && ... args)
		      : // val{std::forward<Arg>(args)...}
		   val{ {static_cast<data_type>(std::forward<Arg>(args))...} }
		    { };

	// Copy ... delete it, this is a base class
	base_ht(const base_ht & t) = delete;
	// ... move
	base_ht(base_ht && t) noexcept =default;
	self_type & operator=(base_ht && t) noexcept = default;
// 	Enable moving a standard array into base_ht class
	base_ht(array_type && t) noexcept: val{std::move(t)}{};
	self_type & operator=(array_type && t) noexcept{
		val=(std::move(t));
		return *this;
	}

	// Is a base class. Virtual destructor ...
	virtual ~base_ht()= default;

	virtual data_type * data() noexcept { return val.data(); }
	virtual const data_type * data() const noexcept { return const_cast<data_type *>(val.data()); }
	// Just to force the interface in case you need to do some calculation or set some timestamps ...

	virtual void operator()() noexcept =0;
	virtual void operator()(const data_type *, std::size_t ) noexcept = 0;
	constexpr std::size_t size() const noexcept { return Size;}

	constexpr data_type & operator[](std::size_t n) noexcept
	      { return val[n]; }

	constexpr const data_type & operator[](std::size_t n) const noexcept
	      { return val[n]; }

    constexpr data_type * begin() noexcept
    { return val.begin(); }

    constexpr const data_type * begin() const noexcept
    { return val.begin(); }

    constexpr data_type * end() noexcept
    { return val.end(); }

    constexpr const data_type * end() const noexcept
    { return val.end(); }

	array_type val;
};


} // namespace message
} // namespace fcutils



#endif // BASE_HT_HPP_ 
