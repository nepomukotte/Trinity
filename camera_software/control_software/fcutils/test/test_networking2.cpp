///
///
/// @brief A new file
///
/// @file test_networking.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 13 ott 2021 
///

#include <iostream>
#include <string>
#include <thread>
#include <boost/program_options.hpp>
#include <networking/client.hpp>
#include <networking/session.hpp>


// Alias: std
namespace ch = std::chrono;
// Alias: boost
namespace ba=boost::asio;
namespace po = boost::program_options;


int main(int ac, char* av[])
{
	// Let's look at the options: ip and port
	// A very basic usage of BOOST program options
	ba::ip::port_type port;
	std::string address;
	// start adding options
	po::options_description desc("test_networking options ");
	desc.add_options()
	        ("help,h", "An help message")
	        ("address,a", po::value<std::string>(&address)->default_value("127.0.0.1"), "the port to connect to")
	        ("port,p", po::value<ba::ip::port_type>(&port)->default_value(20021U), "the port to connect to")
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

	std::cout << " Port: " <<port << '\n';
	std::cout << " Address: " <<address << '\n';

	using tcp_client_type=fcutils::networking::base_client< fcutils::networking::tcp, false, fcutils::networking::session >;
	using udp_client_type=fcutils::networking::base_client< fcutils::networking::udp, false, fcutils::networking::session >;

	try
  {
	boost::asio::io_context io_context;
// Construct a signal set registered for process termination.
	ba::signal_set signals(io_context, SIGINT, SIGTERM);

	// Start an asynchronous wait for one of the signals to occur.
	signals.async_wait([](const boost::system::error_code& ec, int signal_number){
		std::cout << " Signal received: " << signal_number <<'\n';
		if(!ec) {
			std::cout << " Signal handler. No error \n";
		} else {
			std::cout << " signal handler. Error : " << ec.message() << std::endl;
		}
	});
	// A consumer is needed, we must use the one defined in the session
	tcp_client_type::session_type::consumer_type print_data;
    // Create a client on both UDP and TCP/IP
    tcp_client_type ctcp(io_context,print_data);
    udp_client_type cudp(io_context,print_data);
    ba::ip::tcp::endpoint eptcp{boost::asio::ip::make_address(address),port};
    ba::ip::tcp::endpoint eptcp2{boost::asio::ip::make_address(address),static_cast<ba::ip::port_type>(port+1U)};
    ba::ip::udp::endpoint epudp{boost::asio::ip::make_address(address),port};
    ba::ip::udp::endpoint epudp2{boost::asio::ip::make_address(address),static_cast<ba::ip::port_type>(port+1U)};
    ctcp.async_connect(eptcp);
    ctcp.async_connect(eptcp2);
    cudp.async_connect(epudp);
    cudp.async_connect(epudp2);
    ctcp.set_option(eptcp,boost::asio::socket_base::reuse_address{true});
    ctcp.set_option(eptcp2,boost::asio::detail::socket_option::boolean<SOL_SOCKET, SO_REUSEPORT>{true});
    cudp.set_option(epudp,boost::asio::socket_base::reuse_address{true});
    cudp.set_option(epudp2,boost::asio::detail::socket_option::boolean<SOL_SOCKET, SO_REUSEPORT>{true});
    ctcp.start();
    cudp.start();
    // A pointer to the member function run(). Is needed because is an oveloaded function
    ba::io_context::count_type (ba::io_context::*pf_run)(){&ba::io_context::run};
    // Create an launch a thread executing io_context.run()
    std::thread t{pf_run,&io_context};
    std::string s;
    std::cout << "\n Enter a string to broadcast: " << std::endl;
    while(std::cin>>s){
    	ctcp.write(s.c_str(),s.size());
    	ctcp.write(s.c_str(),s.size(),tcp_client_type::endpoint_type(boost::asio::ip::make_address(address),port));
    	cudp.write(s.c_str(),s.size());
    	cudp.write(s.c_str(),s.size(),udp_client_type::endpoint_type(boost::asio::ip::make_address(address),port));
        std::cout << "\n Enter a string to broadcast: " << std::endl;
    }
    ctcp.stop();
    cudp.stop();
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception: " << e.what() << "\n";
  }

  return 0;
}


