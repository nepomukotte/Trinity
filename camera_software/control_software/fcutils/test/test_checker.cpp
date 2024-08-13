///
///
/// @brief A new file
///
/// @file test_checker.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 28 ott 2021 
///

#include <checker/true_checker.hpp>
#include <checker/words_checker.hpp>
#include <iostream>
#include <string>


// Alias: fcutils
namespace fcc=fcutils::checker;

int main(int ac, char* av[])
{

	std::string msg;
	fcc::true_checker<> always_true{'f'};
	fcc::bytes_checker<3> fab_checker{'a',{'f','a','b'}};
//	fcc::bytes_completer<3> fab_completer{'f','a','b'};
	fcc::words_checker<3,std::size_t,fcc::decode_header,char,char> fab_completer{'a',{'f','a','b'},fcc::decode_header<bool,std::size_t>{5,4}};

 	std::cout << " Enter a message to check: \n";
 	while(std::getline(std::cin,msg)){
 		for( auto c: msg ){
 			std::cout << "\n Char: " << c << ", Check: " << std::boolalpha << always_true(&c)<<','<<fab_checker(&c)<<" Completer: " << fab_completer(&c) ;
 		}
 		always_true.reset();
 		fab_checker.reset();
 		std::cout << "\n Message: " << msg << ", Check: " << std::boolalpha << fab_checker(msg.data(),msg.size()) <<" Completer: " << fab_completer(msg.data(),msg.size()) ;
 		fab_checker.reset();
 		fab_completer.reset();
 	 	std::cout << "\n Enter a message to check: \n";
 	}
 	std::cout << std::endl;

  return 0;
}




