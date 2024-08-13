///
///
/// @brief A general purpose header or trailer with a consumer
///
/// @file general_ht.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 10 mar 2022 
///

#ifndef GENERAL_HT_HPP_
#define GENERAL_HT_HPP_

#include <type_traits>

#include "fcutils/message/base_ht.hpp"

namespace fcutils{
namespace message{

enum class GENERAL_HT_OPTION: bool {AT_END=false,AT_BEGIN=true};

// Use SFINAE to chek if the Size is large anough to store a time_t value
template<typename Data, std::size_t Size, typename Consumer, GENERAL_HT_OPTION AtBegin=GENERAL_HT_OPTION::AT_BEGIN,
		typename std::enable_if<sizeof(
				typename std::conditional<
//					std::is_invocable_v< decltype(&Consumer::operator()),Consumer *, const Data *, std::size_t >,
					std::is_invocable_v< Consumer, const Data *, std::size_t >,
//					typename std::invoke_result< decltype(&Consumer::operator()),Consumer *, const Data *, std::size_t>::type ,
					typename std::invoke_result< Consumer, const Data *, std::size_t>::type ,
					typename std::invoke_result< decltype(&Consumer::consume   ), Consumer *, const Data *, std::size_t>::type
				>::type
				) <= Size
				>::type * = nullptr>
struct general_ht: public base_ht<Data,Size> {
	using self_type		=general_ht<Data,Size,Consumer,AtBegin>;
	using base_type		=base_ht<Data,Size>;
	using data_type		=typename base_type::data_type;
	using data_ptr_type	=typename base_type::data_ptr_type;
	using consumer_type	=Consumer;
	using value_type 	=typename std::conditional<
			// Check if the Consumer has an operator(), if so, use the operator() return type, otherwise use the consume() method and return type
			// C++17 is needed .... Is using the CPonsumer pointer, just in case Consumer is a derived class and the consume and operator() are inherited.
//				std::is_invocable<decltype(&consumer_type::operator()),consumer_type *, const data_ptr_type,std::size_t >::value,
				std::is_invocable<consumer_type, const data_type *,std::size_t >::value,
//				typename std::invoke_result<decltype(&consumer_type::operator()),consumer_type *,const data_ptr_type,std::size_t>::type ,
				typename std::invoke_result<consumer_type,const data_type *,std::size_t>::type ,
				typename std::invoke_result<decltype(&consumer_type::consume)   ,consumer_type *,const data_type *,std::size_t>::type
			>::type;

	// The Offset is calculate in bytes. Is important for the "end" case.
	enum: std::size_t {CONSUMER_RESULT_SIZE=sizeof(value_type), CONSUMER_OFFSET=static_cast<bool>(AtBegin)?0:(Size-CONSUMER_RESULT_SIZE)};

	// std::array constructors and assignments are not defined because it follows the rule of an aggregate type. So we would use a little trick ....
	template <typename ... Arg,
			typename std::enable_if<sizeof...(Arg) <= Size>::type * = nullptr>
		   constexpr general_ht (Arg && ... args)
		      : // val{std::forward<Arg>(args)...}
		   base_type{ {static_cast<data_type>(std::forward<Arg>(args))...} }, value_pos_{reinterpret_cast<value_type *>(base_type::data()+CONSUMER_OFFSET)}
		    { };


	bool is_at_end() const { return !static_cast<bool>(AtBegin); }
	bool is_at_begin() const { return static_cast<bool>(AtBegin); }

	std::size_t value_offset() const { return static_cast<std::size_t>(CONSUMER_OFFSET); }

	value_type latest_value() const {
		return *value_pos_;
	}

	void operator()() noexcept { *value_pos_=consumer_(base_type::data(),base_type::size()); }
	void consume() { operator()();}

	void operator()(const data_type * c, std::size_t length) noexcept { *value_pos_=consumer_(c,length);}
	void consume( const data_type * c, std::size_t length) { operator()(c,length);}

	consumer_type & consumer(){ return consumer_;}

private:
	value_type * value_pos_;
	consumer_type consumer_;
};

} // namespace message
} // namesapce fcutils




#endif // GENERAL_HT_HPP_ 
