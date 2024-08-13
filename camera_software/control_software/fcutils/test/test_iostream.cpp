///
///
/// @brief A new file
///
/// @file test_iostream.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 24 mar 2022 
///
#include <iostream>
#include <boost/program_options.hpp>
#include "fcutils/iostream/iostream.hpp"

namespace po=boost::program_options;

int main(int ac, char* av[])
{
	// Let's look at the options: ip and port
	// A very basic usage of BOOST program options
	int port;
	std::string address;
	// start adding options
	po::options_description desc("test_iostream options ");
	desc.add_options()
	        ("help,h", "An help message")
	        ("address,a", po::value<std::string>(&address)->default_value("127.0.0.1"), "the address to connect to")
	        ("port,p", po::value<int>(&port)->default_value(20021), "the port to connect to")
	;
	// Create the option map
	po::variables_map vm;
	// Parse command line for options
	po::store(po::parse_command_line(ac, av, desc), vm);
	// Store user options into the map and the local variables
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 0;
	}

	std::cout << "Address: "<< address <<", Port: " <<port << '\n';

	fcutils::iostream::ostream_logger logger{"test_iostream"};
	// enable log std::cout into a file
	logger.enable(std::cout);
	fcutils::iostream::net_duplicator_manager net_mgr;
	address+=':';
	address+=std::to_string(port);
	address+=" T";

	net_mgr.enable(address.c_str(),std::cout);
	address[address.size()-1]='U';
	net_mgr.enable(address.c_str(),std::cout);
//	fcutils::iostream::ostream_tcp_duplicator tcp_dup(std::cout, address.c_str(), port);
//	fcutils::iostream::ostream_udp_duplicator udp_dup(std::cout, address.c_str(), port);
	int count{0};
	std::string test;
	while(std::getline(std::cin,test)) {
		std::cout <<" Echo :"<< test << std::endl;
		++count;
		// Increment the file counter every 5 entries
		if((count%5)==0){
			logger.disable(std::cout);
			logger.enable(std::cout);
		}
	}

//	tcp_dup.disable(std::cout);
//	udp_dup.disable(std::cout);

//	net_mgr.disable_all(std::cout);

	return 0;
}


