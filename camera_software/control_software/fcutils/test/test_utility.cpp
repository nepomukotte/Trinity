///
///
/// @brief A new file
///
/// @file test_utility.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 3 feb 2022 
///

#include <iostream>
#include <tuple>
#include <variant>
#include <utility>
#include "utility.hpp"

template< typename ... T> struct test{
        test(T... args): a{args ...}{};
        std::tuple<T ...> a;
};

template< typename Ref, typename Val> void assign_tuple_element(Ref &r, Val v){
        r=v;
}

int main(){

	// Test ntuple utilities, operator <<, apply_modifier
        test a{1,2,3};
        std::cout  << a.a << '\n';
        fcutils::apply_modifier(assign_tuple_element<int,int>,a.a,5,6,7);
        std::cout  << a.a << '\n';
        fcutils::apply_modifier([](int & a, const int b){ a=b;},a.a,8,9,10);
        std::cout  << a.a << '\n';

    // Test finally

        // Create a buffer of 100 words
        auto junk = new char[100];
        // create a finally object that will apply the lambda upon distruction
        auto f=fcutils::finally([junk](){ std::cout << " Deleting the buffer ... "<< std::endl; delete junk;});

    // Test Overload
        std::variant<int, std::string, test<int,int,int> > var(42);

        auto ovl=fcutils::overload{
        			[](int i) { std::cout << "int: " << i << '\n'; },
        			[](const std::string& s) { std::cout << "string: " << s << '\n'; },
        			[](const test<int,int,int> & t) { std::cout << "tuple<int,int,int>: " << t.a << '\n';}
	     };

		std::visit(ovl,var);
		var=" This is a variant utility test";
		std::visit(ovl,var);
		var=a;
		std::visit(ovl,var);
		var=" Now I would like to test the print_vars";
		fcutils::print_vals(std::cout, std::get<std::string>(var), a.a, 3.14, 256);

		int ti;
		std::string ts;
		test<int,int,int> tt{4,5,6};
		auto ovl2= fcutils::overload{
					[&ti](int i) mutable { std::cout << "filling int: " << i << '\n'; ti=i;},
					[&ts](const std::string& s) mutable { std::cout << " filling string: " << s << '\n'; ts=s; },
					[&tt](const test<int,int,int> & t) mutable { std::cout << " filling tuple<int,int,int>: " << t.a << '\n'; tt=t;}
		};
		std::visit(ovl2,var);
	// Test integers
		fcutils::int_t<10>::least test;
		fcutils::int_t<64>::least test2;
		fcutils::int_t<17>::least test3;
		fcutils::int_t<6>::least test4;
// Uncomment the next line to get a compiler error
//		fcutils::int_t<65>::least test5;

		std::cout << " Sizes of least integer with 10,64,17 and 6 bits:";
		fcutils::print_vals(std::cout,sizeof(test),sizeof(test2),sizeof(test3),sizeof(test4));
		std::cout << std::endl;
		// A mask of 10 bits, with 8 bit offset: 0x300
		fcutils::mask_t<10,8> mask1;
		// A mask of 64 bits, with 63 bit offset: 0x80000000
		fcutils::mask_t<64,63> mask2;
		// A mask of 3 bits, with 4 bit offset in a 17 bit, at least, integer:0x0070
		fcutils::mask_t<7,4,17> mask3;
		// A mask of 6 bits
		fcutils::mask_t<6> mask4;
		std::cout << " Masks of 10, 64, 17 and 6 bit sizes: "<<std::hex<< std::showbase; // << mask1.mask << std::endl;
		fcutils::print_vals(std::cout,mask1.mask,mask2.mask,mask3.mask,mask4.mask);
		// An or between mask1, a mask of the first byte of an integer of at least 10 bit size and a mask of the second bit of the second byte of at least 10 bit size.
		test = mask1 | fcutils::mask_t<4,0,10>{} | fcutils::mask<10>(6,6);
		std::cout << std::hex << " Mask: "<< test << std::endl;
	// Test extractor
		// An example using CAN J1939
		// Create a CAN J1939 with the PDU (Power Distribution Unit) "default" CAN ID
		fcutils::uint_t<29>::least pduid{0x18efcf50u};
		// Create the mask to extract the source address the PGN and priority
		fcutils::mask_t<8,0,29> maddr;
		fcutils::mask_t<18,0,29> mpgn;
		fcutils::mask_t<3,0,29> mpriority;
		std::cout << " Masks to be used for J1939 extraction: " << std::hex << std::showbase;
		fcutils::print_vals(std::cout,maddr,mpgn,mpriority);
		std::cout << std::endl;
		auto addr=fcutils::extract(pduid,0,maddr.mask);
		auto pgn=fcutils::extract(pduid,8,mpgn.mask);
		auto prio=fcutils::extract(pduid,26,mpriority.mask);
		std::cout << " PDU Address, PGN and priority: " << std::hex << std::showbase;
		fcutils::print_vals(std::cout,addr,pgn,prio);
		std::cout << std::endl;
		// Do the same as above, but using an extractor.
		fcutils::extractor<29,8,0> addrextractor;
		fcutils::extractor<29,18,8> pgnextractor;
		fcutils::extractor<29,3,26> priorityextractor;
		std::cout << " Masks to be used for J1939 extraction (extractor case): " << std::hex << std::showbase;
		fcutils::print_vals(std::cout,addrextractor.get_mask(),pgnextractor.get_mask(),priorityextractor.get_mask());
		std::cout << std::endl;
		addr=addrextractor.extract(pduid);
		pgn=pgnextractor.extract(pduid);
		prio=priorityextractor.extract(pduid);
		std::cout << " PDU Address, PGN and priority (extractor case): " << std::hex << std::showbase;
		fcutils::print_vals(std::cout,addr,pgn,prio);
		std::cout << std::endl;
	// Test inserter
		fcutils::inserter<29,8,0> addrinserter;
		fcutils::inserter<29,18,8> pgninserter;
		fcutils::inserter<29,3,26> priorityinserter;
		pduid=0u;
		addrinserter.insert(pduid,addr);
		pgninserter.insert(&pduid,pgn);
		priorityinserter.insert(pduid,prio);
		std::cout << " CAN ID: " << std::hex << std::showbase << pduid << std::endl;
		// Test std::byte iostream
		std::cout << "Enter a byte (in hex format): ";
		std::byte bt;
		while(std::cin >> std::hex >> bt >> std::dec){
			std::cout << "\n the byte: " <<bt<<std::endl;
			std::cout << "Enter a byte (in hex format): ";
		}

		return 0;
}


