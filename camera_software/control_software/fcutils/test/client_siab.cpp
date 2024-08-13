#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string_view>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <iterator>
#include <cstdint>
#include <cstddef>
#include <bitset>
#include <errno.h>
#include <fcntl.h>
#include <array>
#include <vector>
#include <algorithm>
#include <iomanip>

#include <sys/ioctl.h>
#include <linux/i2c-dev.h>
#include <linux/i2c.h>
#include <i2c/smbus.h>

#include "commands.h"
#include "posix_message_queue.hpp"

using namespace std;
using fcutils::system::posix_message_queue;

int file = 0;
std::string response = "2020";
bool new_cmd = false;

std::vector<uint32_t> HextoInt(const std::string& hex, int msg_size)
{
	std::vector<uint32_t> vmsg;
	for (int i = 0; i < msg_size; i += 2)
	{
		std::string byteString = hex.substr(i, 2);
		vmsg.push_back(stoi(byteString,0,16));
	}
	return vmsg;
}

void Parse_Run_CMD(const std::string& msg, int msg_size, std::string& response)
{
	std::vector<uint32_t> vCMD;

	if ((msg_size%2) > 0){
		cout << "Length of message mismatch. Please enter even number of hex letters." << endl;
		return;
	}
	
	vCMD = HextoInt(msg, msg_size);

	if (vCMD.size() != 8){
		cout << "Length of message mismatch. Please enter 8 bytes or 16 hex letters." << endl;
		return;
	}

	Process_CMD(file, vCMD, response);
	std::cout << "response: " << response << endl;
}

void read_mq(sigval_t sig)
{
	int numRead{-1};
	// Retrive the posix message queue object from the sigval_t infos.
	posix_message_queue* mqp=static_cast<posix_message_queue*>(sig.sival_ptr);
	auto msg_size=mqp->get_msgsize();
	auto msg=new char[msg_size];
	auto mq_name=mqp->get_name();
	mqp->notify(read_mq);
	while ((numRead = mqp->receive(msg, msg_size)) > 0)
	{
		std::cout << " Queue: " << mq_name <<", message received: " << msg << ", of size: " << numRead << std::endl;
		Parse_Run_CMD(msg, numRead, response);
	}
	new_cmd = true;
	delete [] msg;
}

int main()
{
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::cout << "Staring at: " << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << std::endl;

	file = open("/dev/i2c-0", O_RDWR);
	if (file < 0)
	{
		std::cout << "Error in opening file ..." << std::endl;
		exit(1);
	}
	
	std::string mqname="/writetoSIAB";
	std::cout << "Creating and opening the POSIX message queue: " << mqname << " for listening to commands ..." << std::endl;
	posix_message_queue rq(mqname.c_str());
	rq.open();

	mqname="/writetoRC";
	std::cout << "Creating and opening the POSIX message queue: " << mqname << " for responding to master ..." << std::endl;
	posix_message_queue wq(mqname.c_str());
	wq.open();

	rq.notify(read_mq);

	while(true){
		if(new_cmd){
			//wq.send(response);
			//response.clear();
			new_cmd = false;
		}
		usleep(100000);
	}

	return 0;
}