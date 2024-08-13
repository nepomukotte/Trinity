///
///
/// @brief A new file
///
/// @file ansi_format.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 31 gen 2022 
///

#ifndef ANSI_FORMAT_HPP_
#define ANSI_FORMAT_HPP_
#include <iosfwd>
#include "ansi_detail.hpp" // for parameter_type, enum class csi, par2type, operator <<


namespace fcutils{
namespace terminal{
namespace ansi{

enum class format: parameter_type {
	TXT_RESET 		= 0,
	TXT_BOLD		= 1,
	TXT_DIM			= 2,
	TXT_ITALIC		= 3,
	TXT_UNDERLINED	= 4,
	TXT_BLINKSLOW	= 5,
	TXT_BLINKFAST	= 6,
	TXT_REVERSE		= 7,
	TXT_HIDDEN		= 8
};

std::ostream & operator<<(std::ostream & o, const format & f){
	return o << static_cast<trailer_type>(csi::ESCAPE) <<
			static_cast<trailer_type>(csi::SQUAREBRACKET) <<
			static_cast<parameter_type>(f) <<
			static_cast<trailer_type>(csi::SGR);
}


//const par2type<parameter_type, static_cast<parameter_type>(format::TXT_RESET) > TXT_RESET;
//const par2type<parameter_type, static_cast<parameter_type>(format::TXT_BOLD) > TXT_BOLD;
//const par2type<parameter_type, static_cast<parameter_type>(format::TXT_DIM) > TXT_DIM;
//const par2type<parameter_type, static_cast<parameter_type>(format::TXT_ITALIC) > TXT_ITALIC;
//const par2type<parameter_type, static_cast<parameter_type>(format::TXT_UNDERLINED) > TXT_UNDERLINED;
//const par2type<parameter_type, static_cast<parameter_type>(format::TXT_BLINKSLOW) > TXT_BLINKSLOW;
//const par2type<parameter_type, static_cast<parameter_type>(format::TXT_BLINKFAST) > TXT_BLINKFAST;
//const par2type<parameter_type, static_cast<parameter_type>(format::TXT_REVERSE) > TXT_REVERSE;
//const par2type<parameter_type, static_cast<parameter_type>(format::TXT_HIDDEN) > TXT_HIDDEN;
//


} // namespace ansi

} // namesspace terminal
} // namespace fcutils




#endif // ANSI_FORMAT_HPP_ 
