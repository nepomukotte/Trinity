///
///
/// @brief A new file
///
/// @file server.h
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 13 ott 2021 
///

#ifndef SERVER_HPP_
#define SERVER_HPP_

#include <utility> // for std::move
#include <memory> // for std::make_shared
#include <type_traits> // for std::is_same, std::conditional
#include <iostream>
#include <boost/asio.hpp>

namespace fcutils{
namespace networking{

namespace ba=boost::asio;
namespace baip=ba::ip;

using tcp=baip::tcp;
using udp=baip::udp;
using udp_acceptor=ba::basic_socket_acceptor<udp>;

// A type selector for the acceptor_type. From C++11 we can use the std::conditional!!!!!
// Just return the second template type
template <bool IsUdp, class T> struct acceptor_selector{
	using type=typename T::acceptor;
};
// Specialization for UDP, return the
template<> struct acceptor_selector<false,udp>{
	using type=udp_acceptor;
};

// Primary template
template < class Protocol, bool Echo, template<typename,bool> class Session > class base_server{};

// TCP/IP specialization
template <  bool Echo, template<typename,bool> class Session > class base_server<tcp,Echo,Session>
{
public:
	using protocol_type=tcp;
	using session_type=Session<protocol_type,Echo>;
	using acceptor_type=typename protocol_type::acceptor;
	using endpoint_type=typename protocol_type::endpoint;
	using socket_type=typename protocol_type::socket;
	using address_type=ba::ip::address;
	using checker_type=typename session_type::checker_type;
	using consumer_type=typename session_type::consumer_type;

	// Accept connection on a local port
	base_server(ba::io_context& io_context, int port, consumer_type & consumer, const checker_type & checker=checker_type{})
		: acceptor_(io_context, endpoint_type(protocol_type::v4(), port)), consumer_{consumer}, checker_{checker}
	{
		boost::asio::socket_base::reuse_address option;
		acceptor_.get_option(option);
		std::cout << " Reuse_address option: " << std::boolalpha << option.value() << std::endl;

		do_accept();
	}

private:
	void do_accept()
	{
		acceptor_.async_accept(
				[this](boost::system::error_code ec, socket_type socket)
				{
				if (!ec)
				{
					std::make_shared<session_type>(std::move(socket),consumer_,checker_)->start();
				}
				std::cout << " re-running do_accept ..." << std::endl;
				do_accept();
        });
		std::cout << " Async_accept launched ..." << std::endl;
  }

	acceptor_type acceptor_;
	consumer_type & consumer_;
	checker_type checker_;
};

// UDP specialization
template <  bool Echo, template<typename,bool> class Session > class base_server<udp,Echo,Session>
{
public:
	using protocol_type=udp;
	using session_type=Session<protocol_type,Echo>;
	using endpoint_type=typename protocol_type::endpoint;
	using socket_type=typename protocol_type::socket;
	using address_type=ba::ip::address;
	using consumer_type=typename session_type::consumer_type;
	using checker_type= typename session_type::checker_type;

	base_server(ba::io_context& io_context, int port, consumer_type & consumer, const checker_type& checker=checker_type{})
    : io_context_(io_context), port_(port), consumer_{consumer}, checker_{checker}
  {
    do_accept();
  }

  void do_accept()
  {
	  io_context_.post(
		  [this](){
			std::make_shared<session_type>(socket_type(io_context_,endpoint_type(udp::v4(),port_)),consumer_,checker_)->start();
	  	  }
	  );
  }


private:
  ba::io_context & io_context_;
  udp::endpoint endpoint_;
  int port_;
  consumer_type & consumer_;
  checker_type checker_;
};


} // namespace networking
} // namespace fcutils



#endif // SERVER_HPP_ 
