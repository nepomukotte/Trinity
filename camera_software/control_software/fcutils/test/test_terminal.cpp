///
///
/// @brief A new file
///
/// @file test_terminal.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 31 gen 2022 
///

#include <iostream>
#include <string>
#include "terminal/ansi.hpp"

namespace fcta=fcutils::terminal::ansi;

int main(){
	std::cout <<
			fcta::color::FG_BLACK << fcta::color::BG_WHITE << " Background white, foreground black \n" <<
			fcta::color::FG_RED << fcta::color::BG_GREEN << " Background green, foreground red \n" <<
			fcta::color::FG_GREEN << fcta::color::BG_BLUE << " Background blue, foreground green \n" <<
			fcta::color::FG_YELLOW << fcta::color::BG_BLACK << " Background black, foreground yellow \n" <<
			fcta::color::FG_BLUE << fcta::color::BG_MAGENTA << " Background magenta, foreground blue \n" <<
			fcta::color::FG_MAGENTA << fcta::color::BG_YELLOW << " Background yellow, foreground magenta \n" <<
			fcta::color::FG_CYAN << fcta::color::BG_RED << " Background red, foreground cyan \n" <<
			fcta::color::FG_WHITE << fcta::color::BG_CYAN << " Background cyan, foreground white \n" <<
			fcta::color::FG_DEFAULT << fcta::color::BG_DEFAULT << " Background default, foreground default \n" <<
			fcta::format::TXT_BOLD << " Bold \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_DIM << " Dim \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_ITALIC << " Italic \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_UNDERLINED << " Underlined \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_BLINKSLOW << " Blink slow \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_BLINKFAST << " Blink fast \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_REVERSE << " Reverse \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_HIDDEN << " Hidden \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_RESET << " Reset \n";
	std::cout << " Enter return to go home and overwrite ...\n";
	std::string junk;
	std::getline(std::cin,junk);
	std::cout << fcta::CURSORHOME;
	std::cout <<
			fcta::color::FG_DEFAULT << fcta::color::BG_DEFAULT << " Background default, foreground default \n" <<
			fcta::color::FG_WHITE << fcta::color::BG_CYAN << " Background cyan, foreground white \n" <<
			fcta::color::FG_CYAN << fcta::color::BG_RED << " Background red, foreground cyan \n" <<
			fcta::color::FG_MAGENTA << fcta::color::BG_YELLOW << " Background yellow, foreground magenta \n" <<
			fcta::color::FG_BLUE << fcta::color::BG_MAGENTA << " Background magenta, foreground blue \n" <<
			fcta::color::FG_YELLOW << fcta::color::BG_BLACK << " Background black, foreground yellow \n" <<
			fcta::color::FG_GREEN << fcta::color::BG_BLUE << " Background blue, foreground green \n" <<
			fcta::color::FG_RED << fcta::color::BG_GREEN << " Background green, foreground red \n" <<
			fcta::color::FG_BLACK << fcta::color::BG_WHITE << " Background white, foreground black \n" <<
			fcta::format::TXT_BOLD << " Bold \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_DIM << " Dim \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_ITALIC << " Italic \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_UNDERLINED << " Underlined \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_BLINKSLOW << " Blink slow \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_BLINKFAST << " Blink fast \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_REVERSE << " Reverse \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_HIDDEN << " Hidden \n" << fcta::format::TXT_RESET <<
			fcta::format::TXT_RESET << " Reset \n";
	int n,m;
	std::cout << " Enter the number of column to move the cursor: ";
	std::cin >> n;
	std::getline(std::cin,junk);
	fcta::CURSORATCOLUMN(std::cout,n);
	std::cout << fcta::color::FG_RED << fcta::color::BG_GREEN << " Background green, foreground red " << fcta::color::FG_DEFAULT << fcta::color::BG_DEFAULT << '\n' ;
	std::cout << " Enter the number of column and line to move the cursor: ";
	std::cin >> n >> m;
	std::getline(std::cin,junk);
	fcta::CURSORXYPOSITION(std::cout,m,n);
	std::cout << fcta::color::FG_YELLOW << fcta::color::BG_BLUE << fcta::format::TXT_BOLD << " HERE I AM !!! " << fcta::format::TXT_RESET << fcta::color::FG_DEFAULT << fcta::color::BG_DEFAULT << '\n' ;
	std::cout << " Enter return to clear display ...\n";
	std::cin.clear();
	std::getline(std::cin,junk);
	fcta::clear();
	std::cout << " Enter the number of column and line to move the cursor: ";
	std::cin >> n >> m;
	std::getline(std::cin,junk);
	fcta::gotoxy(n,m);
	std::cout << " Enter the number of column and line to reset the terminal: ";
	std::cin >> n >> m;
	std::getline(std::cin,junk);
	fcta::resize(n,m);
	return 0;
}


