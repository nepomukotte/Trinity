///
///
/// @brief A new file
///
/// @file test_status_msg.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 6 ott 2021 
///

#include <boost/asio.hpp>
#include <boost/program_options.hpp>
#include <checker/words_checker.hpp>
#include <networking/client.hpp>
#include <networking/server.hpp>
#include <networking/session.hpp>
#include <iostream>
#include <string>
#include <chrono>
#include <ctime>
#include <cstring>
#include <thread>
#include <cstddef>


// Alias: std
namespace ch = std::chrono;
// Alias: boost
namespace ba=boost::asio;
namespace po = boost::program_options;

//// User defined literal operator to convert integer into char
//constexpr char operator"" _int_to_char(){
//
//}
namespace fcc=fcutils::checker;
namespace fcn=fcutils::networking;

const std::size_t GCCHEADER_SIZE{8},GCCSTATUS_MAX_SIZE{104};
using data_type=char;
using key_type=char;
using checker_return_type=std::size_t;
const bool DEBUG{false}; // Set this to true for more DEBUG printout from the checker
template<typename Data> using gcc_consumer_type=fcn::print_byte<Data>;
template<typename Data> using gcc_checker_type=fcc::words_checker<GCCHEADER_SIZE,checker_return_type,fcc::decode_header,Data,key_type,DEBUG>;
template<typename Protocol, bool Echo> using gcc_session_type=fcn::session<Protocol, Echo, gcc_consumer_type, gcc_checker_type >;
using gcc_client_type=fcn::base_client< fcn::tcp, false, gcc_session_type >;

int main(int ac, char* av[])
{

	// Let's look at the options: ip and port
	// A very basic usage of BOOST program options
	std::string address;
	ba::ip::port_type port;
	int delay{0};
	bool timestamp{false};
	// start adding options
    po::options_description desc("test_status_msg options ");
    desc.add_options()
        ("help,h", "An help message")
        ("ip", po::value<std::string>(&address)->default_value("192.168.1.1"), "the IP address to talk to")
        ("port,p", po::value<ba::ip::port_type>(&port)->default_value(20021U), "the port to connect to")
		("delay,d", po::value<int>(&delay)->default_value(30), "the delay between broadcast (seconds)")
		("timestamp,t", po::bool_switch(&timestamp))
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

    std::cout << " IP address: " << address <<":"<<port<<", delay (sec): "<<delay<< '\n';
    std::string s=(timestamp)?"timestamp":"message counter";
    std::cout << "Using "<<s<<"as status message \n";

    // Here the checker for words of size char, storable in a checker map with a char key.
	// Here the number of words to be checked is the size of the header
    // This object will be copied into a client or a server.
	gcc_checker_type<data_type> gcc_checker{
		'g',
		{ // The header byte constant values to be checked for only 7 of 8 bytes will be used,
			static_cast<data_type>(0xAA), // The header ID, four bytes 0xAA5555AA
			static_cast<data_type>(0x55),
			static_cast<data_type>(0x55),
			static_cast<data_type>(0xAA),
			static_cast<data_type>(0x00), // The message type: 0x0, for received msg; 0x1 for status, i.e. transmitted, message
			static_cast<data_type>(0x00), // Mandatory empty bytes, two bytes
			static_cast<data_type>(0x00)
		},
		fcc::decode_header<bool,std::size_t>{8,8} // The number of words to be checked upon extraction of the size from the offset.
	}; // end of words_checker constructor

	// We need a consumer.
	// The consumer is supposed to be passed by reference.
	// This is because a consumer can be also a container that could be filled with new consumers at run time.
	// For this example we simply use the standard print consumer, used as default in the session parameter list.
	gcc_client_type::session_type::consumer_type print_data;

    // Try to deal with the TCP stream objects
    try{
    	// start with an ASIO context
    	ba::io_context ioc;
    	// Construct a signal set registered for process termination.
    	ba::signal_set signals(ioc, SIGINT, SIGTERM);
    	// Start an asynchronous wait for one of the signals to occur.
    	signals.async_wait([](const boost::system::error_code& ec, int signal_number){
    		std::cout << " Signal received: " << signal_number <<'\n';
    		if(!ec) {
    			std::cout << " Signal handler. No error \n";
    		} else {
    			std::cout << " signal handler. Error : " << ec.message() << std::endl;
    		}
    	});
    	// The client TCP to be used with a session checking for the GCC header
    	gcc_client_type gcc_client{ioc,print_data,gcc_checker};
    	// Start the client
    	gcc_client.start();
    	// Add the address and port and try to connect
    	gcc_client.async_connect(boost::asio::ip::make_address(address),port);
    	// Start the context on another thread
        // A pointer to the member function run(). Is needed because is an oveloaded function
        ba::io_context::count_type (ba::io_context::*pf_run)(){&ba::io_context::run};
        // Create an launch a thread executing io_context.run()
        std::thread t{pf_run,&ioc};
       	// create a time point converted to time_t (C style) for an easier text conversion
       	std::time_t now=ch::system_clock::to_time_t(ch::system_clock::now());
       	// create a time duration for the wait state
       	auto wd=ch::seconds{delay};
       	int count{0};
       	// Prepare a standard header needed by the GCC message
       	const std::size_t GCCHEADER_SIZE{8},GCCSTATUS_MAX_SIZE{104};
       	unsigned char gccheader[GCCHEADER_SIZE]={static_cast<unsigned char>(0xAA),
       			static_cast<unsigned char>(0x55),static_cast<unsigned char>(0x55),
				static_cast<unsigned char>(0xAA),static_cast<unsigned char>(0x01),
				static_cast<unsigned char>(0x00),static_cast<unsigned char>(0x00) };
       	// Just the null termination char
       	char gcctrailer='\0';
       	// create the buffer sequence
       	std::array<ba::const_buffer,2> status_msg;
       	// A pointer to the message
       	char *cnow{nullptr};
       	// The message length
       	std::size_t length{0};
       	while(1){
       		++count;
       		std::cout<< "Sending status message no. "<<count<<" ...\n";
       		if(timestamp){
       			// this is in principle a static array of char ...
       			auto now=ch::system_clock::to_time_t(ch::system_clock::now());
       			cnow=std::ctime(&now);
       			// Remember, std::strlen returns the lenght without the null termination
       			length=cnow?std::strlen(cnow)+1:0;
       			// Set maximum length to 146 bytes
       			length=(length<=GCCSTATUS_MAX_SIZE)?length:GCCSTATUS_MAX_SIZE;
       		} else {
       			cnow=reinterpret_cast<char *>(&count);
       			length=sizeof(count);
       		}
       		// Store the message length into the gcc header
       		gccheader[7]=static_cast<char>(length);
//       		stream<<std::ctime(&now);
//       		stream << count;
       		std::cout << " Writing the message :"<< cnow << ", of length: "<< length <<'\n';
       		// Create the buffer sequence
       		status_msg[0]=ba::buffer(gccheader,GCCHEADER_SIZE);
       		status_msg[1]=ba::buffer(cnow,length);
       		gcc_client.write(status_msg);
           	// Put this thread on wait for 30 seconds
           	std::this_thread::sleep_for(wd);
       	}


    }
    // Catch ASIO exceptions here
    catch (std::exception &e)
        {
            std::cerr << "Exception: " << e.what() << "\n";
        }//

	return 0;
}

