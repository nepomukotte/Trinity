///
///
/// @brief To test fcutils/system utilities
///
/// @file test_system.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 21 feb 2022 
///

#include <iostream>
#include <boost/asio.hpp>
#include "fcutils/system/getip.hpp"
#include "fcutils/networking/udpbeacon.hpp"

int main(){
	std::string dev,msg;
	fcutils::networking::udpbeacon beacon{1s};
	fcutils::networking::udpjsonbeacon jbeacon{1500ms};
	beacon.start();
	jbeacon.start();
	std::cout << " Enter the device name: ";
	while (std::cin>>dev){
		msg=dev;
		msg+=" <" + fcutils::system::getip(dev) + '>';
		msg+=", netmask <" + fcutils::system::getip(dev,fcutils::system::netmask) + '>';
		msg+=", broadcast <" + fcutils::system::getip(dev,fcutils::system::broadcast) + '>';
		msg+='\n';
		std::cout << " Ip addresses for device: " << msg << std::endl;
		beacon.set_message(msg);
		std::cout << " Enter the device name: ";
	}

	beacon.stop();
	jbeacon.stop();
	return 0;
}


