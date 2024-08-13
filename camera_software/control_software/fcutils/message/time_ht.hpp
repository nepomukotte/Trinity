///
///
/// @brief An header or trailer with a time stamp
///
/// @file time_ht.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 9 mar 2022 
///

#ifndef TIME_HT_HPP_
#define TIME_HT_HPP_
#include <ctime>
#include "fcutils/message/base_ht.hpp"

namespace fcutils{
namespace message{

enum class TIME_HT_OPTION: bool {AT_END=false,AT_BEGIN=true};

// Use SFINAE to check if the Size is large anough to stor a time_t value
template<typename Data, std::size_t Size, TIME_HT_OPTION AtBegin=TIME_HT_OPTION::AT_BEGIN, typename std::enable_if<sizeof(std::time_t) <= Size>::type * = nullptr>

struct time_ht: public base_ht<Data,Size> {
	// The Offset is calculate in bytes. Is important for the "end" case.

	enum: std::size_t {CTIME_SIZE=sizeof(std::time_t), TIME_OFFSET=static_cast<bool>(AtBegin)?0:(Size-CTIME_SIZE)};

	using self_type=time_ht<Data,Size,AtBegin>;
	using base_type=base_ht<Data,Size>;
	using data_type=typename base_type::data_type;
	using data_ptr_type=typename base_type::data_ptr_type;

	// std::array constructors and assignments are not defined because it follows the rule of an aggregate type. So we would use a little trick ....
	template <typename ... Arg,
			typename std::enable_if<sizeof...(Arg) <= Size>::type * = nullptr>
		   constexpr time_ht (Arg && ... args)
		      : // val{std::forward<Arg>(args)...}
		   base_type{ {static_cast<data_type>(std::forward<Arg>(args))...} }, time_pos_{reinterpret_cast<std::time_t *>(base_type::data()+TIME_OFFSET)}
		    { };


	bool at_end() const { return !static_cast<bool>(AtBegin); }
	bool at_begin() const { return static_cast<bool>(AtBegin); }

	std::size_t time_offset() const { return static_cast<std::size_t>(TIME_OFFSET); }

	std::time_t latest_timestamp() const {
		return *time_pos_;
	}

	void generate_timestamp() const {
		*time_pos_=std::time(NULL);
	}

	operator std::time_t() const noexcept { return *time_pos_;}

	void operator()() noexcept { generate_timestamp(); }
	void consume() { operator()();}

	void operator()(const data_type * c, std::size_t length) noexcept { generate_timestamp();}
	void consume( const data_type * c, std::size_t length) { operator()(c,length);}

private:
	std::time_t * time_pos_;

};

} // namespace message
} // namespace futils



#endif // TIME_HT_HPP_ 
