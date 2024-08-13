//
// test_string.cpp
//
//  Created on: 19 dic 2020
//      Author: cafagna
//

#include <iostream>
#include <array>
#include <iterator>
#include <cstdint> // for int8_t etc. etc.
#include <cstddef> // for std::byte
#include <bitset>
#include <string_view>

#include "fcutils/utility/utility.hpp" // for operator<< on std::tuple and print vals
#include "fcutils/string/string.hpp"

struct test{
	// We need to be aligned to the largest type ....
	uint8_t i8[2];
	uint16_t i16[3];
};

int main(){

	std::cout<< " Enter a prefix (Integer value):";
	std::string s;
	std::array<char,32> buffer;
	fcutils::string::program_name_and_date pnd;

	while(std::cin>>s){
		// Test hash function
		std::cout << " String: " << s << ", hash: 0x" << std::hex << fcutils::string::generic_hash(s.c_str(), fcutils::string::hash_xor_op) << ", 0x" << fcutils::string::generic_hash(s.c_str(), fcutils::string::hash_sum_op) << std::dec<< std::endl;
		std::cout << " String: " << s << ", hash: 0x" << std::hex << fcutils::string::hash(s.c_str()) << std::dec<< std::endl;
		// Test add_timestamp
		std::cout << '\n' << fcutils::string::add_timestamp(s) << ", " << fcutils::string::add_timestamp<fcutils::string::pos::BACK>(s)<< '\n';
		std::cout << '\n' << fcutils::string::add_timestamp<fcutils::string::pos::FRONT, fcutils::string::conv::CTIMEFMT>(s) << ", "
				<< fcutils::string::add_timestamp<fcutils::string::pos::BACK, fcutils::string::conv::CTIMEFMT>(s)<< '\n';
		// Test program_name_and_date
		std::cout << '\n' << pnd << '\n';
		// Create a string with the ASCI value of the characters composing the string s
		auto temp=fcutils::string::char_to_hex(s.c_str(),s.size());
		std::cout << " Prefix character ASCII value in hex: " << temp << '\n';
		std::cout << " Using an output stream iterator: ";
		std::ostream_iterator<char> oo3{std::cout," "}; // output iterator for filestream
		fcutils::string::char_to_hex(oo3,s.c_str(),s.size());
		std::cout << '\n';
		std::cout << " Convert back the ASCII value in hex: " << fcutils::string::hex_to_char(temp.c_str(),temp.size());
		std::cout << std::endl;
		// Convert the string into an integer
		auto val=std::stoi(s);
		// Create a string with the representation of the integer number in the prefix
		temp=fcutils::string::int_to_hex(val);
		std::cout << " Prefix character representation in hex: " << temp << '\n';
		// Reset the string and fill it with the representation
		auto last=fcutils::string::int_to_hex(buffer.data(),val);
		// Terminate the buffer to print the right number of char ...
		*last=0;
		std::cout << " Saving to an existing buffer: " << buffer.data() << '\n';
		// Because the target buffer is a template, you can use an iterator ... also a stream one
		std::cout << " Using an output stream iterator: ";
		std::ostream_iterator<char> oo{std::cout," "}; // output iterator for filestream
		fcutils::string::int_to_hex(oo,val);
		std::cout << std::endl;
		std::cout << " Using a string: ";
		s+=' ';
		fcutils::string::int_to_hex(s,val);
		std::cout << s << std::endl;
		val=0;
		fcutils::string::hex_to_int(temp.c_str(),val);
		std::cout << " Converting back an integer value: " << val << '\n';
		test my_test{0xfaU,0xb0U,0xFAb1U,0xbac1U,0xeabbu};
		unsigned char *pt{reinterpret_cast<unsigned char*>(&my_test)},*ptend{pt+sizeof(test)};
		std::string my_test_str{"0x"};
		for(auto i=pt;i!=ptend;++i) fcutils::string::int_to_hex(my_test_str,*i,"");
		std::cout << " Structure serialization: " << my_test_str << '\n';
		test my_test_copy{};
		pt=reinterpret_cast<unsigned char*>(&my_test_copy);
		fcutils::string::hex_to_int(pt,my_test_str.c_str(),my_test_str.size());
		std::cout << " Original struct: " << std::hex;
		fcutils::print_vals(std::cout, static_cast<uint16_t>(my_test.i8[0]),static_cast<uint16_t>(my_test.i8[1]), my_test.i16[0], my_test.i16[1], my_test.i16[2]);
		std::cout << " Copied struct: " << std::hex;
		fcutils::print_vals(std::cout, static_cast<uint16_t>(my_test_copy.i8[0]),static_cast<uint16_t>(my_test.i8[1]), my_test_copy.i16[0], my_test_copy.i16[1], my_test_copy.i16[2]);
		std::cout<< "\n Enter a prefix (Integer value):";
	}
	// Test case conversion
	std::string tc{"This Is Mixed"};
	// convert to uppercase and lowercase
	std::cout << "\n To upper case: " << fcutils::string::to_upper(tc) << "\n  To lower case: " << fcutils::string::to_lower(std::move(tc)) << "\n Original : " << tc << std::endl;
	// The same with an array of char
	char tcc[]="This Is Mixed";
	// convert to uppercase and lowercase
	std::cout << " To upper case: " << fcutils::string::to_upper(tcc) << "\n  To lower case: " << fcutils::string::to_lower(tcc) << "\n Original : " << tcc << std::endl;


	return 0;
}

