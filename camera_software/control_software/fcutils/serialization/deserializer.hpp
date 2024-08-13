///
///
/// @brief A simple deserializer
///
/// @file deserializer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 8 feb 2022 
///

#ifndef DESERIALIZER_HPP_
#define DESERIALIZER_HPP_

#include <string>
#include <algorithm>
#include <type_traits>
#include "string/hex_to_int.hpp" // for int_to_hex functions

namespace fcutils{
namespace serialization{


struct deserialize_representation{
	template<typename ToBeDeserialized, typename Data, typename SizeType > auto operator()(ToBeDeserialized *s,
			const Data * source, SizeType length){
		auto pt=reinterpret_cast<unsigned char*>(s);
		return fcutils::string::hex_to_int(pt,source,length);
	}
};

struct just_copy{

	template<typename ToBeDeserialized, typename Data> auto operator()(ToBeDeserialized *s, const Data * source, fcutils::string::type_size<ToBeDeserialized,2ul> length){
		return operator()(s,source,static_cast<std::size_t>(length));
	}

	template<typename ToBeDeserialized, typename Data> auto operator()(ToBeDeserialized *s, const Data * source, std::size_t length){
		// Special case for string
		if constexpr( std::is_same<ToBeDeserialized,std::string>::value && std::is_same<Data,char>::value){
			return s->append(source, length);
		} else {
			return std::copy_n(source,length,s);
		}
	}
};

template <typename Policy> class deserializer{
public:
	using policy_type=Policy;

	template<typename ToBeDeserialized, typename Data>
	auto deserialize(ToBeDeserialized * s, const Data * source, std::size_t length ){
		return deserialize_(s, source, length);
	}

	template<typename ToBeDeserialized, typename Data, std::size_t ByteCharSize=2ul >
	auto deserialize(ToBeDeserialized * s, const Data * source) {
		return deserialize_(s, source, fcutils::string::type_size<ToBeDeserialized,ByteCharSize>{});
	}
private:
	policy_type deserialize_;
};

using deserialize_repr=deserializer<deserialize_representation>;
using deserialize_string=deserializer<just_copy>;

} // namespace serialization
} // namespace fcutils






#endif // DESERIALIZER_HPP_ 
