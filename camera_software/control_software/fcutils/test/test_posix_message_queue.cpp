///
///
/// @brief A new file
///
/// @file test_posix_message_queue.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 14 giu 2021 
///

#include <string>
#include <iostream>

#include "posix_message_queue.hpp"

using namespace std;
using fcutils::system::posix_message_queue;

const std::string MQCS2RC("/mqueue_cs2rc"),MQRC2CS("/mqueue_rc2cs");

void read_mq(sigval_t sig){
	int numRead{-1};
	// Retrive the posix message queue object from the sigval_t infos.
	posix_message_queue* mqp=static_cast<posix_message_queue*>(sig.sival_ptr);
	auto msg_size=mqp->get_msgsize();
	auto msg=new char[msg_size];
	auto mq_name=mqp->get_name();
	mqp->notify(read_mq);
	while ((numRead = mqp->receive(msg, msg_size)) > 0){
		std::cout << " Queue: " << mq_name <<", message received: " << msg << ", of size: " << numRead<< std::endl;
	}
	delete [] msg;
}

int main(){
	std::cout << " Enter the name of the message queue to listen to [/mqueue_cs2rc]: " << std::endl;
	std::string mqname;
	std::getline(std::cin, mqname);
	mqname=mqname.empty()?MQCS2RC:mqname;
	std::cout << " Creating and opening the POSIX message queue with name: " << mqname << std::endl;
	posix_message_queue rq(mqname.c_str());
	rq.open();
	std::cout << " Enter the name of the message queue to write into [/mqueue_rc2cs]: " << std::endl;
	std::getline(std::cin, mqname);
	mqname=mqname.empty()?MQRC2CS:mqname;
	std::cout << " Creating and opening the POSIX message queue with name: " << mqname << std::endl;
	posix_message_queue wq(mqname.c_str());
	wq.open();
	rq.notify(read_mq);
	std::string msg;
	std::cout << "Enter the message to be sent to the Run Control \n";
	while(std::getline(std::cin,msg)){
		std::cout << " Sending message: " << msg <<", of size: "<< msg.size() << std::endl;
		wq.send(msg);
	}
	return 0;
}



