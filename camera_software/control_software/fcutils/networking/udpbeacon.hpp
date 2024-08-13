///
///
/// @brief An udp beacon ... using the old way ...
///
/// @file udpbecaon.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 21 feb 2022 
///

#ifndef UDPBEACON_HPP_
#define UDPBEACON_HPP_

#include <string>
#include <iostream>
#include <iomanip> // for std::quote
#include <cstring> // for std::memset
#include <thread> // for std::thread
#include <shared_mutex> // for std::shared_mutex
#include <chrono>
#include <atomic> // for std::atomic
#include <utility> // for std::forward
#include <netinet/in.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <errno.h> // for program_invocation_short_name

#include "fcutils/concurrency/guarded_thread.hpp"
#include "fcutils/system/getip.hpp"

namespace{
	using namespace std::chrono_literals;
}

namespace fcutils{
namespace networking{

// Message creator, to be passed to the udpbeacon
struct standard_xml_msg{
	int port{9000};
	std::string operator()(const std::string & device="eth0"){
		std::string msg;
			msg="<UDPbeacon>\n<APP>";
			// Uset the application invocation name
			msg+=program_invocation_short_name;
			msg+="</APP>";
			msg+="<URI>tcp:";
			// Get the IP address
			msg+=fcutils::system::getip(device);
			msg+=':';
			msg+=std::to_string(port);
			msg+="</URI>";
			msg+="</UDPbeacon>";
			return msg;
	}
};

// Message creator, to be passed to the udpbeacon
struct standard_json_msg{
	int port{9000};
	std::string operator()(const std::string & device="eth0"){
		std::string msg;
		msg="{\"UDPbeacon\": {\n\"APP\": \"";
		// Uset the application invocation name
		msg+=program_invocation_short_name;
		msg+="\",\"URI\": \"";
		msg+="tcp:";
		// Get the IP address
		msg+=fcutils::system::getip(device);
		msg+=':';
		msg+=std::to_string(port);
		msg+="\"}}";
		return msg;
	}
};


template <typename Duration=std::chrono::seconds, typename DefMsgGenerator=fcutils::networking::standard_xml_msg>
class base_udpbeacon{
public:
	using socket_type=	int;
	using port_type=	int;
	using duration_type=Duration;
	using msg_generator_type=DefMsgGenerator;
	using message_type=	std::string;
	using self_type=	base_udpbeacon<Duration,DefMsgGenerator>;

	// An guarded_thread will take care of joining the thread upon destruction
	using thread_type=	fcutils::concurrency::guarded_thread;

	template<typename D=duration_type> base_udpbeacon(D && d, const std::string & dev="eth0",port_type port=8888):
			wait_{std::forward<duration_type>(std::chrono::duration_cast<duration_type>(d))},
			device_{dev},
			port_{port} {};

	~base_udpbeacon(){
		// Stop the beacon. It is a stopping operation ...
		do_stop();
		// Close the socket if any
		if(sock_>0) ::close(sock_);
	}

	void start()
	{
		if(started_.load()) return;
//		std::cout << " Starting the UDP beacon ... " << std::endl;
		started_.store(true);
		// create the socket
		create_socket();
		// create the standard message if msg is not set
		if(msg_.empty()) create_msg();
		// launch the broadcast
		b_thread_=std::thread(&self_type::do_broadcast,this);
	}

	void stop()
	{
		if(!started_.load()) return;
//		std::cout << " Stopping the UDP beacon ... " << std::endl;
		stop_thread_=std::thread(&self_type::do_stop,this);
	}

	template <typename M=message_type> void set_message(M && m){
	    std::shared_lock lock(mutex_);
		msg_=std::forward<message_type>(m);
	}

	template <typename D=duration_type> void set_wait(D && d){
	    std::shared_lock lock(mutex_);
		wait_=std::forward<duration_type>(d);
	}

private:
	void create_msg(){
		msg_generator_type g;
		msg_=g(device_);
	}

	void create_socket(){
    // Open a socket
		sock_ = ::socket(PF_INET, SOCK_DGRAM, IPPROTO_UDP);
		if (sock_<=0) {
			std::cout<< "UDP beacon: Error opening the socket" << std::endl;
			return;
		}
		// Set socket options
		// Enable broadcast
		int enable=1;
		int ret=::setsockopt(sock_, SOL_SOCKET, SO_BROADCAST, &enable, sizeof(enable));
		if (ret) {
			std::cout<< "UDP beacon: Error setting the broadcast option" << std::endl;
			return;
		}
		// Enable reuse address and port
		ret=::setsockopt(sock_, SOL_SOCKET, SO_REUSEADDR, &enable, sizeof(enable));
		if (ret) {
			std::cout<< "UDP beacon: Error setting the reuse address option" << std::endl;
			return;
		}
		ret=::setsockopt(sock_, SOL_SOCKET, SO_REUSEPORT, &enable, sizeof(enable));
		if (ret) {
			std::cout<< "UDP beacon: Error setting the reuse port option" << std::endl;
			return;
		}
		// Configure the port and ip we want to send to
	    sockaddr_in addr; // Make an endpoint
	    std::memset(&addr, 0, sizeof(addr));
	    // set the family
	    addr.sin_family = AF_INET;
	    // Set the broadcast IP address and port
	    ipaddr_=fcutils::system::getip(device_,fcutils::system::broadcast);
	    ::inet_pton(AF_INET,ipaddr_.c_str() , &addr.sin_addr);
	    addr.sin_port = ::htons(port_);
	    ::connect(sock_,reinterpret_cast<sockaddr*>(&addr), sizeof(addr));
	}

	void do_stop(){
		if(!started_.load()) return;
		// Set the started flag to false
		started_.store(false);
		// waits for the broadcast thread
		b_thread_.join();
		// Thats's all
	}

	void do_broadcast(){
		while(started_.load()){
			// Just wait for the wating time and broadcast
			std::this_thread::sleep_for(wait_);
			broadcast();
		}
	}

	void broadcast (){
	    // lock the message and broadcast it
	    std::shared_lock lock(mutex_);
	    // Broadcast the message
	    int ret = ::send(sock_, msg_.c_str(), msg_.size(), 0);
	    if (ret<0) {
	        std::cout<< "beacon: Error sending the message" << std::endl;
	        return;
	    }
	}

	std::shared_mutex mutex_;
	std::atomic_bool started_{false};
	port_type port_;
	socket_type sock_{-1};
	message_type msg_;
	thread_type b_thread_;
	thread_type stop_thread_;
	duration_type wait_{std::chrono::duration_cast<duration_type>(5s)}; // 5 seconds is the default duration
	std::string ipaddr_;
	std::string device_;
};

using udpbeacon=base_udpbeacon<>;
using udpjsonbeacon=base_udpbeacon<std::chrono::seconds, fcutils::networking::standard_json_msg>;

} // namespace networking
} // namespace fcutils



#endif // UDPBEACON_HPP_ 
