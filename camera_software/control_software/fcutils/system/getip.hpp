///
///
/// @brief Retrive the local IP number, either IP4 and IP6
///
/// @file getip.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 21 feb 2022 
///

#ifndef GETIP_HPP_
#define GETIP_HPP_

#include <cstdio>
#include <string>
#include <iostream>
#include <type_traits> // for std::is_same
#include <ifaddrs.h>
#include <netdb.h>
#include <sys/socket.h>
#include "fcutils/utility/finally.hpp"

namespace fcutils{
namespace system{

// Define types to select different type of address extraction. Namely: address, for the IP address; netmask, for the netmask; broadcastaddres, for the boradcast address.
class get_address {};
class get_netmask {};
class get_broadcastaddress {};
// The actual objects for tag dispatch
inline constexpr get_address address;
inline constexpr get_netmask netmask;
inline constexpr get_broadcastaddress broadcast;

template<typename GetType=get_address > auto getip( const std::string & dev="eth0", GetType gt=fcutils::system::address, int family=AF_INET) {
	// Start with an X ip address ...
	std::string temp{" no device "};
    ifaddrs *ifaddr{nullptr};
    int s;
    const std::size_t IPADDRSIZE{25};
    char ipaddr[IPADDRSIZE];

    if (::getifaddrs(&ifaddr) == -1) {
    	std::cout << " getip: Error in getifaddr " << std::endl;
    	return std::string{};
    }
    // free the ifaddr once done
    auto f=fcutils::finally([ifaddr](){ if(ifaddr) ::freeifaddrs(ifaddr);});
    // The pointer to the actual socket address to be queried
    sockaddr *addr=nullptr;
    // Walk through linked list, maintaining head pointer so we can free list later
    for (auto ifa = ifaddr; ifa != nullptr; ifa = ifa->ifa_next) {
    	// Check what kind of address is requested. If the standard address, the netmask or the broadcast address
    	if constexpr(std::is_same<GetType,get_address>::value) {
    		addr=ifa->ifa_addr;
    	} else {
    		if constexpr(std::is_same<GetType,get_netmask>::value) {
    			addr=ifa->ifa_netmask;
    		} else {
    			addr=ifa->ifa_broadaddr;
    		}
    	}
        if (addr){
        	// Check if the device name is the same
        	if(dev==ifa->ifa_name){
        		// For an AF_INET* interface address, get the address
        		if (family == ifa->ifa_addr->sa_family) {
        			s = ::getnameinfo(addr,
        					(family == AF_INET) ? sizeof(sockaddr_in) : sizeof(sockaddr_in6),
        					ipaddr, IPADDRSIZE,
							nullptr, 0, NI_NUMERICHOST);
        			if (s != 0) {
        				std::cout << " getip: Error in getnameinfo " << std::endl;
        				return std::string{};
        			}
        			temp=ipaddr;
        		} // if (family ...
        	} // if(dev ...
        } // if(ifa->ifa_addr)

    } // for (ifa = ifaddr, n = 0; ifa != nullptr; ifa = ifa->ifa_next, ++n)

//    if (ifaddr!=NULL) ::freeifaddrs(ifaddr);
    return temp;
}

} // namespace system
} // namespace fcutils


#endif // GETIP_HPP_ 
