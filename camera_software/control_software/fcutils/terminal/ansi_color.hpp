///
///
/// @brief An enum class with the ANSI color definition
///
/// @file color.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 31 gen 2022 
///

#ifndef ANSI_COLOR_HPP_
#define ANSI_COLOR_HPP_
#include <iosfwd>
#include "ansi_detail.hpp" // for parameter_type, enum class csi, par2type, operator <<

namespace fcutils{
namespace terminal{
namespace ansi{

enum class color: parameter_type {
	FG_BLACK	= 30,
	FG_RED      = 31,
	FG_GREEN    = 32,
	FG_YELLOW	= 33,
	FG_BLUE     = 34,
	FG_MAGENTA	= 35,
	FG_CYAN		= 36,
	FG_WHITE	= 37,
	FG_DEFAULT  = 39,
	BG_BLACK	= 40,
	BG_RED      = 41,
	BG_GREEN    = 42,
	BG_YELLOW	= 43,
	BG_BLUE     = 44,
	BG_MAGENTA	= 45,
	BG_CYAN		= 46,
	BG_WHITE	= 47,
	BG_DEFAULT  = 49,
};

std::ostream & operator<<(std::ostream & o, const color & c){
	// insure that the output format is dec ...
	std::ios_base::fmtflags ff{o.flags()};
	// Set the format to decimal
	o.setf (std::ios_base::dec , std::ios::basefield);
	o << static_cast<trailer_type>(csi::ESCAPE) <<
			static_cast<trailer_type>(csi::SQUAREBRACKET) <<
			static_cast<parameter_type>(c) <<
			static_cast<trailer_type>(csi::SGR);
	o.flags(ff);
	return o;
}

//const par2type<parameter_type, static_cast<parameter_type>(color::FG_BLACK) > FG_BLACK;
//const par2type<parameter_type, static_cast<parameter_type>(color::FG_RED) > FG_RED;
//const par2type<parameter_type, static_cast<parameter_type>(color::FG_GREEN) > FG_GREEN;
//const par2type<parameter_type, static_cast<parameter_type>(color::FG_YELLOW) > FG_YELLOW;
//const par2type<parameter_type, static_cast<parameter_type>(color::FG_BLUE) > FG_BLUE;
//const par2type<parameter_type, static_cast<parameter_type>(color::FG_MAGENTA) > FG_MAGENTA;
//const par2type<parameter_type, static_cast<parameter_type>(color::FG_CYAN) > FG_CYAN;
//const par2type<parameter_type, static_cast<parameter_type>(color::FG_WHITE) > FG_WHITE;
//const par2type<parameter_type, static_cast<parameter_type>(color::FG_DEFAULT) > FG_DEFAULT;
//const par2type<parameter_type, static_cast<parameter_type>(color::BG_BLACK) > BG_BLACK;
//const par2type<parameter_type, static_cast<parameter_type>(color::BG_RED) > BG_RED;
//const par2type<parameter_type, static_cast<parameter_type>(color::BG_GREEN) > BG_GREEN;
//const par2type<parameter_type, static_cast<parameter_type>(color::BG_YELLOW) > BG_YELLOW;
//const par2type<parameter_type, static_cast<parameter_type>(color::BG_BLUE) > BG_BLUE;
//const par2type<parameter_type, static_cast<parameter_type>(color::BG_MAGENTA) > BG_MAGENTA;
//const par2type<parameter_type, static_cast<parameter_type>(color::BG_CYAN) > BG_CYAN;
//const par2type<parameter_type, static_cast<parameter_type>(color::BG_WHITE) > BG_WHITE;
//const par2type<parameter_type, static_cast<parameter_type>(color::BG_DEFAULT) > BG_DEFAULT;



} // namespace ANSI
} // namespace terminal

} // namespace fcutils




#endif // ANSI_COLOR_HPP_ 
