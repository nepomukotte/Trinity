///
///
/// @brief A new file
///
/// @file test_serialport.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 27 ott 2021 
///
#include <iostream>
#include <string>
#include <algorithm>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/asio.hpp>
#include <boost/optional.hpp>
#include <concurrency/guarded_thread.hpp>

#include <serialport/talk_to_serial.hpp>

// Alias: boost
namespace ba=boost::asio;
namespace po = boost::program_options;

int main(int ac, char* av[])
{
	// Let's look at the options: ip and port
	// A very basic usage of BOOST program options
	boost::optional<std::string> port;
	// start adding options
	po::options_description desc("test_serialport options ");
	desc.add_options()
	        ("help,h", "An help message")
	        ("port,p", po::value<boost::optional<std::string>>(&port), "the port to connect to")
	;
	// Create the option map
	po::variables_map vm;
	// Parse command line for options
	po::store(po::parse_command_line(ac, av, desc), vm);
	// Store user options into the map and the local variables
	po::notify(vm);

	if (vm.count("help")||!port) {
		std::cout << " To test the program you can create a dummy serial port using: socat \n socat -d -d pty,raw,echo=0 pty,raw,echo=0 \n"
				<< " Then run test_serialport using the ports created by socat, typically: /dev/pts/1 and /dev/pts/2. \n"
				<< " So in one terminal enter: ./test_serialport -p /dev/pts/1\n"
				<< " So in one other terminal enter: ./test_serialport -p /dev/pts/2\n\n";
		std::cout << desc << "\n";
		return 0;
	}

	std::cout << " Using serial port: " <<*port << '\n';

	try
  {

    boost::asio::io_context io_context;
   	// Create the listening process on the serial port
    // Use the
    auto ttsp=fcutils::serialport::talk_to_serial<>::get(io_context,*port);
    ttsp->start_read();
    // Launch the io_context on a separate thread
     // A pointer to the member function run()
 	ba::io_context::count_type (ba::io_context::*pf_run)(){&ba::io_context::run};
 	// Create an launch a thread executing io_context.run()
 	fcutils::concurrency::guarded_thread read_thread{pf_run,&io_context};
 	std::string msg;
 	std::cout << " Enter a message to write into the serial port: \n";
 	while(std::getline(std::cin,msg)){
 		// Add a carriage return
 		msg+='\n';
 		// remeber to write also the null termination ... (we need to add an extra byte)
 		ttsp->write(msg.c_str(),msg.size()+1);
 	 	std::cout << " Enter a message to write into the serial port: \n";
 	}
 	io_context.stop();
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception: " << e.what() << "\n";
  }

  return 0;
}



