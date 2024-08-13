///
///
/// @brief A new file
///
/// @file tcp_connection_tester.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 24 mar 2022 
///

#ifndef TCP_CONNECTION_TESTER_HPP_
#define TCP_CONNECTION_TESTER_HPP_
#include <boost/asio.hpp>
#include <chrono>
#include <atomic>
#include <iosfwd>

namespace fcutils{
namespace networking{

class tcp_connection_tester{
public:
	tcp_connection_tester(const int t=100): wtime_{100}, started_{false} {}

	bool test_connection(const char *ipaddress, const int port)  {
		boost::asio::io_service ios;
		boost::asio::ip::tcp::endpoint ep(boost::asio::ip::address::from_string(ipaddress), port);
		boost::asio::ip::tcp::socket sock(ios);
		started_=false;
		sock.async_connect(ep,
				[this](boost::system::error_code ec){ connect_handler(ec);});
		boost::asio::steady_timer t{ios};
		t.expires_after(wtime_);
		t.async_wait(
				[this](boost::system::error_code ec){ deadline_handler(ec);});
		ios.run();
		sock.close();
		return started_.load();
	}

	// In this case the seocket and the endpoint are external and the socket is not closed after the test.
	template < typename Socket, typename Endpoint > bool test_connection(Socket & sock , Endpoint & ep)  {
		started_=false;
		sock.async_connect(ep,
				[this](boost::system::error_code ec){ connect_handler(ec);});
		boost::asio::steady_timer t{sock.get_executor()};
		t.expires_after(wtime_);
		t.async_wait(
				[this](boost::system::error_code ec){ deadline_handler(ec);});
		sock.get_executor().run();
//		sock.close();
		return started_.load();
	}

	void connect_handler(const boost::system::error_code& error)
	{
	  if (!error)
	  {
//		  std::cout << "Socket is opened." << std::endl;
		  started_=true;
		  // Connect succeeded.
	  } else {
//		  std::cout << "Error in open socket: " << error << std::endl;
		  started_=false;
	  }
	}


	void deadline_handler(const boost::system::error_code &e){
		switch(started_.load()) {
		case true:
	//		std::cout << "Connect done " << std::endl;
			break;
		case false:
	//		std::cout << "Connect timeout " << e << std::endl;
			break;
		}
	}

private:
	std::atomic<bool> started_;
	std::chrono::milliseconds wtime_;
};

} // networking
} // fcutils



#endif // TCP_CONNECTION_TESTER_HPP_ 
