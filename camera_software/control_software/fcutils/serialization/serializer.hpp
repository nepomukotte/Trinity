///
///
/// @brief A simple serializer.
///
/// @file serializer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 8 feb 2022 
///

#ifndef SERIALIZER_HPP_
#define SERIALIZER_HPP_
#include <string>
#include "string/int_to_hex.hpp" // for int_to_hex functions

namespace fcutils{
namespace serialization{


struct serialize_representation{
	template<typename ToBeSerialized> void operator()(std::string & t, const ToBeSerialized *s, std::size_t length=sizeof(ToBeSerialized)){
		auto *pt{reinterpret_cast<const unsigned char*>(s)},*ptend{pt+length};
		for(auto i=pt;i!=ptend;++i) fcutils::string::int_to_hex(t,*i,"");
	}
};

struct just_append{
	template<typename ToBeSerialized> void operator()(std::string & t,const ToBeSerialized *s, std::size_t length=sizeof(ToBeSerialized)){
		t.append(s, length);
	}
};

template <typename Policy> class serializer{
public:
	using data_type=std::string;
	using policy_type=Policy;

	void clear(){ data_.clear(); }
	void reserve( std::string::size_type s) { data_.reserve(s); }
	template<typename ToBeSerialized > auto serialize(const ToBeSerialized * s, std::size_t length=sizeof(ToBeSerialized), const std::string & prefix=""){
		data_+=prefix;
		serialize_(data_,s,length);
	}
	// Of course this works if type T is concatenable to a std::string ...
	template<typename T> auto & operator+=( const T & t){ data_+=t; return *this; }
	auto const & data() const { return data_;}
	auto size() { return data_.size(); }
	auto const * c_str() const { return data_.c_str(); }
private:
	data_type data_;
	policy_type serialize_;
};

using serialize_repr=serializer<serialize_representation>;
using serialize_string=serializer<just_append>;

} // namespace serialization
} // namespace fcutils



#endif // SERIALIZER_HPP_ 
