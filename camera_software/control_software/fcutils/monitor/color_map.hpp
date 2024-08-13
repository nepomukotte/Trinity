///
///
/// @brief A new file
///
/// @file color_map.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 26 mar 2022 
///

#ifndef COLOR_MAP_HPP_
#define COLOR_MAP_HPP_
#include <map>
#include <iostream>
#include <algorithm>

#include "fcutils/terminal/ansi.hpp"
#include "fcutils/monitor/xmlconst.hpp"
#include "fcutils/string/hash.hpp"
#include "fcutils/string/change_case.hpp"

namespace fcutils{
namespace monitor{

enum class fg_color: std::size_t {
	FG_RED=fcutils::string::hash(xml::slowcontrol::variable::value::SVRED.data()),
	FG_GREEN=fcutils::string::hash(xml::slowcontrol::variable::value::SVGREEN.data()),
	FG_BLUE=fcutils::string::hash(xml::slowcontrol::variable::value::SVBLUE.data())
};

enum class bg_color: std::size_t {
	BG_DEFAULT=fcutils::string::hash(xml::slowcontrol::variable::value::SVNONE.data()),
	BG_RED=fcutils::string::hash(xml::slowcontrol::variable::value::SVRED.data()),
	BG_GREEN=fcutils::string::hash(xml::slowcontrol::variable::value::SVGREEN.data()),
	BG_BLUE=fcutils::string::hash(xml::slowcontrol::variable::value::SVBLUE.data())
};

using string_to_color_type = std::map<std::size_t,fcutils::terminal::ansi::color>;

const string_to_color_type STRING_TO_COLOR {
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVRED.data()),fcutils::terminal::ansi::color::FG_RED},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVGREEN.data()),fcutils::terminal::ansi::color::FG_GREEN},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVBLUE.data()),fcutils::terminal::ansi::color::FG_BLUE},
//																		(xml::slowcontrol::variable::value::SBLUE,BLUE)
//																		(xml::slowcontrol::variable::value::SBLUE,BLACK)
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVDEFAULTKEY.data()),fcutils::terminal::ansi::color::FG_BLUE }
};
string_to_color_type STRING_TO_BG_COLOR{
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVNONE.data()),fcutils::terminal::ansi::color::BG_DEFAULT},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVRED.data()),fcutils::terminal::ansi::color::BG_RED},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVGREEN.data()),fcutils::terminal::ansi::color::BG_GREEN},
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVBLUE.data()),fcutils::terminal::ansi::color::BG_BLUE},
//																		(xml::slowcontrol::variable::value::SBLUE,BLUE)
//																		(xml::slowcontrol::variable::value::SBLUE,BLACK)
	{fcutils::string::hash(xml::slowcontrol::variable::value::SVDEFAULTKEY.data()),fcutils::terminal::ansi::color::BG_DEFAULT}
};

} // namespace monitor
} // namespace fcutils



#endif // COLOR_MAP_HPP_ 
