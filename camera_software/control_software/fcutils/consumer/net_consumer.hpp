///
///
/// @brief A network consumer
///
/// @file tcp_udp_consumer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 13 mar 2022 
///

#ifndef NET_CONSUMER_HPP_
#define NET_CONSUMER_HPP_
#include <string>
#include <iosfwd>

#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ip/udp.hpp>
#include <boost/asio/io_service.hpp>

#include "fcutils/consumer/consumer_base.hpp"
#include "fcutils/checker/true_checker.hpp"

namespace fcutils{
namespace consumer{

template <typename Data=unsigned char, template <typename D> class Checker = fcutils::checker::true_checker, typename Protocol=boost::asio::ip::tcp>
class net_consumer: public consumer_base<Data,Checker> {
public:
	using self_type			=net_consumer<Data,Checker,Protocol>;
	using protocol_type		=Protocol;
	using consumer_base_type=consumer_base<Data,Checker>;
	using checker_type 		=typename consumer_base_type::checker_type;
	using data_type			=typename consumer_base_type::data_type;
	using endpoint_type		=typename protocol_type::endpoint;
	using socket_type		=typename protocol_type::socket;
	using port_type			=boost::asio::ip::port_type;

	template<typename C=checker_type> net_consumer(const char * ip="127.0.0.1",port_type port=514, C && checker=fcutils::checker::true_checker<data_type>{}):
			consumer_base_type{std::forward<checker_type>(checker)}, ip_address_{ip}, port_{port},
			ep_{boost::asio::ip::address::from_string(ip_address_.c_str()), port_}{};

	self_type & operator()(const data_type * c, std::size_t length) override {
		// Create the socket, open it and bind it to the local endpoint
		socket_type sock(ios_);
		boost::system::error_code ec;
		sock.connect(ep_,ec);
		if(!ec){
			// Socket connected
			sock.send( boost::asio::buffer(c, length));
			consumer_base_type::update_stats(length);
		} else {
			std::cout << " net_consumer: Socket not connected  ... "<< ec.message() << ", " << ep_<< std::endl;
		}
		sock.close();
		return *this;
	}

private:
	std::string ip_address_;
	port_type port_;
	boost::asio::io_service ios_;
	endpoint_type ep_;
};

template<typename Data=unsigned char, template <typename D> class Checker=fcutils::checker::true_checker>
using tcp_consumer = net_consumer<Data,Checker,boost::asio::ip::tcp>;
template<typename Data=unsigned char, template <typename D> class Checker=fcutils::checker::true_checker>
using udp_consumer = net_consumer<Data,Checker,boost::asio::ip::udp>;


} // namespace consumer
} // namespace fcutils




#endif // NET_CONSUMER_HPP_ 
