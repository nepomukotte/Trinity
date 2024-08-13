///
///
/// @brief A new file
///
/// @file test_networking.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 13 ott 2021 
///

#include <iostream>
#include <boost/program_options.hpp>
#include "fcutils/networking/server.hpp"
#include "fcutils/networking/session.hpp"
#include "fcutils/networking/tcp_connection_tester.hpp"

// Alias: std
namespace ch = std::chrono;
// Alias: boost
namespace ba=boost::asio;
namespace po = boost::program_options;

class udp_server
{
public:
  udp_server(ba::io_context& io_context, short port)
    : socket_(io_context, ba::ip::udp::endpoint(ba::ip::udp::v4(), port))
  {
    do_receive();
  }

  void do_receive()
  {
    socket_.async_receive_from(
        boost::asio::buffer(data_, max_length), sender_endpoint_,
        [this](boost::system::error_code ec, std::size_t bytes_recvd)
        {
          if (!ec && bytes_recvd > 0)
          {
            do_send(bytes_recvd);
          }
          else
          {
            do_receive();
          }
        });
  }

  void do_send(std::size_t length)
  {
	  data_[0]='F';
	  socket_.async_send_to(
        boost::asio::buffer(data_, length), sender_endpoint_,
        [this](boost::system::error_code /*ec*/, std::size_t /*bytes_sent*/)
        {
          do_receive();
        });
  }

private:
  ba::ip::udp::socket socket_;
  ba::ip::udp::endpoint sender_endpoint_;
  enum { max_length = 1024 };
  char data_[max_length];
};


int main(int ac, char* av[])
{
	// Let's look at the options: ip and port
	// A very basic usage of BOOST program options
	int port;
	// start adding options
	po::options_description desc("test_networking options ");
	desc.add_options()
	        ("help,h", "An help message")
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

	std::cout << " Port: " <<port << '\n';

	using tcp_server_type=fcutils::networking::base_server< fcutils::networking::tcp, false, fcutils::networking::session >;
	using udp_server_type=fcutils::networking::base_server< fcutils::networking::udp, false, fcutils::networking::session >;
//	using udp_server_type=udp_server;

	fcutils::networking::tcp_connection_tester tcp_test;
	std::cout << std::boolalpha << " Testing the tcp connection ... " << tcp_test.test_connection("127.0.0.1",port) << std::endl;

	try
  {

    boost::asio::io_context io_context;
    // Create a consumer
    tcp_server_type::session_type::consumer_type print_data;
   	// Create the listening servers on both UDP and TCP/IP
    tcp_server_type stcp(io_context, port, print_data);
    udp_server_type sudp(io_context, port, print_data);

    io_context.run();
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception: " << e.what() << "\n";
  }

  return 0;
}


