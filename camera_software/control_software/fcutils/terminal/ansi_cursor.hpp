///
///
/// @brief A new file
///
/// @file cursor_movement.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 31 gen 2022 
///

#ifndef ANSI_CURSOR_HPP_
#define ANSI_CURSOR_HPP_
#include <iosfwd>
#include "ansi_detail.hpp" // trailer_type, enum class csi

namespace fcutils{
namespace terminal{
namespace ansi{

enum class cursor: trailer_type {
	CURSORUP='A',
	CURSORDOWN='B',
	CURSORRIGHT='C',
	CURSORLEFT='D',
	CURSORNEXTLINE='E',
	CURSORPREVIOUSLINE='F',
	CURSORATCOLUMN='G',
	CURSORXYPOSITION='f',
	CURSORONELINEUP='M',
	CURSORSAVEPOSITION='s',
	CURSORRESTOREPORITION='u',
	CURSORHOME='H',
};

// Stream of the cursor class
std::ostream & operator<<(std::ostream & o, const cursor & c){
	return o << static_cast<trailer_type>(c);
}

//template<trailer_type Command> struct command{ };

//const command<static_cast<trailer_type>(cursor::CURSORHOME)> CURSORHOME;
const par2type<trailer_type,static_cast<trailer_type>(cursor::CURSORHOME),false> CURSORHOME;
const par2type<trailer_type,static_cast<trailer_type>(cursor::CURSORONELINEUP),false> CURSORONELINEUP;

escape_sequence<cursor,cursor::CURSORUP> CURSORUP;
escape_sequence<cursor,cursor::CURSORDOWN> CURSORDOWN;
escape_sequence<cursor,cursor::CURSORRIGHT> CURSORRIGHT;
escape_sequence<cursor,cursor::CURSORLEFT> CURSORLEFT;
escape_sequence<cursor,cursor::CURSORNEXTLINE> CURSORNEXTLINE;
escape_sequence<cursor,cursor::CURSORPREVIOUSLINE> CURSORPREVIOUSLINE;
escape_sequence<cursor,cursor::CURSORATCOLUMN> CURSORATCOLUMN;
escape_sequence<cursor,cursor::CURSORXYPOSITION> CURSORXYPOSITION;

void up(std::size_t line, std::ostream & o= std::cout) {CURSORUP(o,line); }
void down(std::size_t line, std::ostream & o= std::cout) {CURSORDOWN(o,line); }
void forward( std::size_t column, std::ostream & o= std::cout) {CURSORRIGHT(o,column); }
void back( std::size_t column, std::ostream & o= std::cout) { CURSORLEFT(o,column); }
void gotoline(std::size_t line, std::ostream & o = std::cout) { fcutils::terminal::ansi::CURSORNEXTLINE(o, line); }
void backtoline(std::size_t line, std::ostream & o = std::cout) { fcutils::terminal::ansi::CURSORPREVIOUSLINE(o, line); }
void gotox(std::size_t column, std::ostream & o = std::cout) { fcutils::terminal::ansi::CURSORATCOLUMN(o, column); }
void gotoxy(std::size_t column, std::size_t line, std::ostream & o = std::cout){ fcutils::terminal::ansi::CURSORXYPOSITION(o,column,line);}

} // namespace ansi
} // namesapce cursor

} // namesapce fcutils




#endif // ANSI_CURSOR_HPP_ 
