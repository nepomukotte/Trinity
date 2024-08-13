///
///
/// @brief A new file
///
/// @file outformat_map.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 5 apr 2022 
///

#ifndef OUTFORMAT_MAP_HPP_
#define OUTFORMAT_MAP_HPP_
#include <ios>
#include <map>

#include "fcutils/monitor/xmlconst.hpp"
#include "fcutils/string/hash.hpp"
#include "fcutils/string/change_case.hpp"

namespace fcutils{
namespace monitor{

using string_to_outformat_type = std::map<std::size_t,std::ios_base::fmtflags>;
string_to_outformat_type STRING_TO_OUTFORMAT{
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVDEC.data()),std::ios_base::dec},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVHEX.data()),std::ios_base::hex},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVDEFAULTKEY.data()),std::ios_base::dec}
};

} // namespace monitor
} // namespace fcutils



#endif // OUTFORMAT_MAP_HPP_ 
