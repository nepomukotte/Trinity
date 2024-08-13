///
///
/// @brief A new file
///
/// @file inouttype_map.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 5 apr 2022 
///

#ifndef INOUTTYPE_MAP_HPP_
#define INOUTTYPE_MAP_HPP_
#include <map>
#include <variant>

#include "fcutils/monitor/xmlconst.hpp"
#include "fcutils/string/change_case.hpp"
#include "fcutils/string/hash.hpp"
#include "fcutils/type_traits/type_list.hpp"

namespace fcutils{
namespace monitor{

/// The output types. Please note that these are also the index for the value_type. So correspondence must be preserved
enum class inouttype: std::size_t	{
	INT,FLOAT,DOUBLE,BOOL//,CHAR,STRING
};
/// The typelist corresponsing to the enum class inouttype
using typelist_type = fcutils::type_traits::type_list<int,float,double,bool>;

/// The actual value type. Correspondence between the variant declaration and the outtype must be preserved to use the outtype as variant index
using value_type = std::variant<int,float,double,bool>;

using string_to_inouttype_type = std::map<std::size_t,inouttype>;
// Create a global object to map the key to the enum
string_to_inouttype_type STRING_TO_INOUTTYPE{
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVINT.data()),inouttype::INT},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVFLOAT.data()),inouttype::FLOAT},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVDOUBLE.data()),inouttype::DOUBLE},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVBOOL.data()),inouttype::BOOL},
//	{xml::slowcontrol::variable::value::SBOOL,inouttype::CHAR},
//	{xml::slowcontrol::variable::value::SBOOL,inouttpye::STRING},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVDEFAULTKEY.data()),inouttype::DOUBLE}
};

// Helper functions to generate variant or enable alternatives in the variant according to the iouttype enumerator
// Just return a std::variant inizialized to the type indexed by the inouttype value
auto enable_alternative(inouttype t){
	switch(t){
	case inouttype::INT:
		return value_type{typelist_type::type<static_cast<std::size_t>(inouttype::INT)>{}};
		break;
	case inouttype::FLOAT:
		return value_type{typelist_type::type<static_cast<std::size_t>(inouttype::INT)>{}};
		break;
	case inouttype::DOUBLE:
		return value_type{typelist_type::type<static_cast<std::size_t>(inouttype::INT)>{}};
		break;
	case inouttype::BOOL:
		return value_type{typelist_type::type<static_cast<std::size_t>(inouttype::INT)>{}};
		break;
	default:
		return value_type{typelist_type::type<static_cast<std::size_t>(inouttype::INT)>{}};
		break;
	};
}

template <typename Variant> auto enable_alternative(inouttype t, Variant & v){
	switch(t){
	case inouttype::INT:
		return v= typelist_type::type< static_cast<std::size_t>(inouttype::INT) >{};
		break;
	case inouttype::FLOAT:
		return v= typelist_type::type< static_cast<std::size_t>(inouttype::FLOAT) >{};
		break;
	case inouttype::DOUBLE:
		return v= typelist_type::type< static_cast<std::size_t>(inouttype::DOUBLE) >{};
		break;
	case inouttype::BOOL:
		return v= typelist_type::type< static_cast<std::size_t>(inouttype::BOOL) >{};
		break;
	default:
		return v= typelist_type::type< static_cast<std::size_t>(inouttype::INT) >{};
		break;
	};
}

template<typename Val> auto set_alternative(inouttype t, Val val){
	switch(t){
	case inouttype::INT:
		return value_type{static_cast<typelist_type::type<static_cast<std::size_t>(inouttype::INT)>>(val)};
		break;
	case inouttype::FLOAT:
		return value_type{static_cast<typelist_type::type<static_cast<std::size_t>(inouttype::FLOAT)>>(val)};
		break;
	case inouttype::DOUBLE:
		return value_type{static_cast<typelist_type::type<static_cast<std::size_t>(inouttype::DOUBLE)>>(val)};
		break;
	case inouttype::BOOL:
		return value_type{static_cast<typelist_type::type<static_cast<std::size_t>(inouttype::BOOL)>>(val)};
		break;
	default:
		return value_type{static_cast<typelist_type::type<static_cast<std::size_t>(inouttype::INT)>>(val)};
		break;
	};
}

template <typename Variant, typename Val> auto set_alternative(inouttype t, Variant & v, Val val){
	switch(t){
	case inouttype::INT:
		return v= static_cast< typelist_type::type< static_cast<std::size_t>(inouttype::INT) > >(val);
		break;
	case inouttype::FLOAT:
		return v= static_cast< typelist_type::type< static_cast<std::size_t>(inouttype::FLOAT) > >(val);
		break;
	case inouttype::DOUBLE:
		return v= static_cast< typelist_type::type< static_cast<std::size_t>(inouttype::DOUBLE) > >(val);
		break;
	case inouttype::BOOL:
		return v= static_cast< typelist_type::type< static_cast<std::size_t>(inouttype::BOOL) > >(val);
		break;
	default:
		return v= static_cast< typelist_type::type< static_cast<std::size_t>(inouttype::INT) > >(val);
		break;
	};
}

} // namespace monitor
} // namespace fcutils


#endif // INOUTTYPE_MAP_HPP_ 
