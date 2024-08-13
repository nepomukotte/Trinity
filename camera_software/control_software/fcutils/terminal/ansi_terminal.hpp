///
///
/// @brief A new file
///
/// @file ansi_terminal.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 31 gen 2022 
///

#ifndef ANSI_TERMINAL_HPP_
#define ANSI_TERMINAL_HPP_
#include "ansi_detail.hpp" // for parameter_type, enum class csi, par2type, operator <<


namespace fcutils{
namespace terminal{
namespace ansi{

enum class terminal: trailer_type {
	ERASEDISPLAY='J',
	ERASELINE='K',
	SCROLLUP='S',
	SCROLLDOWN='T',
	RESIZETRAILER='t',
	RESIZEHEADER='8'
};

// Stream of the terminal class
std::ostream & operator<<(std::ostream & o, const terminal & c){
	return o << static_cast<trailer_type>(c);
}

const par2type<trailer_type,static_cast<trailer_type>(terminal::ERASEDISPLAY),false> ERASEDISPLAY;
escape_sequence<terminal,terminal::RESIZETRAILER> RESIZEDISPLAY;


void clear(std::ostream & o=std::cout) {o << fcutils::terminal::ansi::ERASEDISPLAY;}
void resize(std::size_t height, std::size_t width, std::ostream & o = std::cout){
	fcutils::terminal::ansi::RESIZEDISPLAY(o,static_cast<trailer_type>(fcutils::terminal::ansi::terminal::RESIZEHEADER),height,width);
}

void resize_clear(std::size_t height, std::size_t width, std::ostream & o = std::cout){
	fcutils::terminal::ansi::resize(height,width,o);
	fcutils::terminal::ansi::clear();
}

} // namespace ansi
} // namespace terminal
} // namesapce fcutils




#endif // ANSI_TERMINAL_HPP_ 
