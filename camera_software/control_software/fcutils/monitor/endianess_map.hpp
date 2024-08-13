///
///
/// @brief A new file
///
/// @file endianess_map.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 5 apr 2022 
///

#ifndef ENDIANESS_MAP_HPP_
#define ENDIANESS_MAP_HPP_
#include <map>

#include "fcutils/monitor/xmlconst.hpp"
#include "fcutils/string/change_case.hpp"
#include "fcutils/string/hash.hpp"

namespace fcutils{
namespace monitor{

/// The endianess types
enum class endianess {BIGENDIAN,LITTLEENDIAN};
using string_to_endianess_type = std::map<std::size_t, endianess>;
string_to_endianess_type STRING_TO_ENDIANESS{
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVBIGENDIAN.data()),endianess::BIGENDIAN},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVLITTLEENDIAN.data()),endianess::LITTLEENDIAN},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVDEFAULTKEY.data()),endianess::BIGENDIAN}
};

} // namespace monitor
} // namespace fcutils



#endif // ENDIANESS_MAP_HPP_ 
