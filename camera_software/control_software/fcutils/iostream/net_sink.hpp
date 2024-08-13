///
///
/// @brief A network sink using the BOOST iostream
///
/// @file net_sink.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 24 mar 2022 
///

#ifndef NET_SINK_HPP_
#define NET_SINK_HPP_

#include <iostream>
#include <string>
#include <type_traits>
#include <boost/iostreams/tee.hpp> // tee_stream & tee_device
#include <boost/iostreams/stream.hpp> // iostream
#include <boost/asio.hpp> // UDP & ASIO facilities

#include "fcutils/networking/tcp_connection_tester.hpp"

namespace fcutils{
namespace iostream{

// General template
template <class Transport> class net_sink {};

// Specialized for udp transport
template <> class net_sink<boost::asio::ip::udp>: public boost::iostreams::sink {
public:
	typedef boost::asio::ip::udp transport_type;
	net_sink(const char * ip="127.0.0.1",int port=514): ip_address_(ip), port_(port){};
	std::streamsize write(const char *s, std::streamsize n)
	{
		boost::asio::io_service ios;
		boost::asio::ip::udp::endpoint ep(boost::asio::ip::address::from_string(ip_address_.c_str()), port_);
		boost::asio::ip::udp::socket sock(ios);
		sock.open(boost::asio::ip::udp::v4());
		sock.send_to( boost::asio::buffer(s, n),ep);
		sock.close();
		return n;
	}
	const char* get_address() const {return ip_address_.c_str();}
	int get_port() const {return port_;}
	// In the UDP case the test of the connection is always true.
	bool test_connection() { return true; }
	enum: bool { is_udp=std::is_same_v<transport_type,boost::asio::ip::udp>, is_tcp=std::is_same_v<transport_type,boost::asio::ip::tcp> };

private:
	std::string ip_address_;
	int port_;
};


template <> class net_sink<boost::asio::ip::tcp>: public boost::iostreams::sink {
public:
	typedef boost::asio::ip::tcp transport_type;
	net_sink(const char * ip="127.0.0.1",int port=514): ip_address_(ip), port_(port){};
	std::streamsize write(const char *s, std::streamsize n)
	{
		boost::asio::io_service ios;
		boost::asio::ip::tcp::endpoint ep(boost::asio::ip::address::from_string(ip_address_.c_str()), port_);
		boost::asio::ip::tcp::socket sock(ios);
		boost::system::error_code ec;
		std::streamsize nsent{0};
		sock.connect(ep,ec);
		if (ec)
		{
			// Silently do nothing ....
//			std::cout << " TCP net sink, error connecting the socked on address: "<< ip_address_<<", port: " << port_ << ", error: "<< ec.message() << '\n';

		} else {
			nsent = sock.send( boost::asio::buffer(s, n));
		}
		sock.close();
		// it must be the same as the write to avoid troubles ...
		return n;
	}
	const char* get_address() const {return ip_address_.c_str();}
	int get_port() const {return port_;}
	bool test_connection() { fcutils::networking::tcp_connection_tester t; return t.test_connection(ip_address_.c_str(), port_); }
	enum: bool { is_udp=std::is_same_v<transport_type,boost::asio::ip::udp>, is_tcp=std::is_same_v<transport_type,boost::asio::ip::tcp> };
private:
	std::string ip_address_;
	int port_;
};

typedef net_sink<boost::asio::ip::udp> udp_sink;
typedef net_sink<boost::asio::ip::tcp> tcp_sink;


} // namespace iostream
} // namesapce fcutils




#endif // NET_SINK_HPP_ 
