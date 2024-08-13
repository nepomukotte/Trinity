#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <unistd.h>
#include <iomanip>

#include "fcutils/utility/utility.hpp"
#include "fcutils/string/string.hpp"
#include "posix_message_queue.hpp"

#define CMD_SETUP		0x01
#define CMD_TURN_ON		0x02
#define CMD_TURN_OFF	0x03
#define CMD_GET_STATE	0x04
#define CMD_GET_USAGE	0x05

using namespace std;
using fcutils::system::posix_message_queue;

int file = 0;
std::string response = "6060";
std::string PDU_DIR = "/home/cherenkov/Programs/control_software/pdu";
bool new_cmd = false;

void Can0Setup(std::string &response)
{
	system("sudo ip link set can0 type can bitrate 250000 && sudo ip link set up can0 &");
}

void TurnOnChannel(int ChannelNumber, std::string &response)
{
	std::string command = "cd " + PDU_DIR + " && python3 pdu.py can0 ch" + std::to_string(ChannelNumber) + " on &";
	system(command.c_str());
}

void TurnOffChannel(int ChannelNumber, std::string &response)
{
	std::string command = "cd " + PDU_DIR + " && python3 pdu.py can0 ch" + std::to_string(ChannelNumber) + " off &";
	system(command.c_str());
}

void GetChannelState(int ChannelNumber, std::string &response)
{
	std::string command = "cd " + PDU_DIR + " && python3 pdu.py can0 ch" + std::to_string(ChannelNumber) + " state &";
	system(command.c_str());
}

void GetChannelUsage(int ChannelNumber, std::string &response)
{
	std::string command = "cd " + PDU_DIR + " && python3 pdu.py can0 ch" + std::to_string(ChannelNumber) + " usage &";
	system(command.c_str());
}

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

void Process_CMD(std::vector<uint32_t>& vCMD, std::string &response)
{
	uint32_t CMD_ID  = vCMD.at(0);
	uint32_t Channel = vCMD.at(1);

	switch (CMD_ID)
	{
		case CMD_SETUP:
			Can0Setup(response);
			break;
		case CMD_TURN_ON:
			TurnOnChannel(Channel, response);
			break;
		case CMD_TURN_OFF:
			TurnOffChannel(Channel, response);
			break;
		case CMD_GET_STATE:
			GetChannelState(Channel, response);
			break;
		case CMD_GET_USAGE:
			GetChannelUsage(Channel, response);
			break;
		default:
			break;
	}
	response="6060";
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

	std::cout << "message: " << msg << endl;
	Process_CMD(vCMD, response);
	std::cout << "response: " << response << endl;
}

void read_mq(sigval_t sig)
{
	int numRead{-1};
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
    std::cout << "Staring PDU at: " << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << std::endl;

	std::string mqname="/writetoPDU";
	std::cout << "Creating and opening the POSIX message queue: " << mqname << " for listening to commands ..." << std::endl;
	posix_message_queue rq(mqname.c_str());
	rq.open();

	mqname="/writetoRC";
	std::cout << "Creating and opening the POSIX message queue: " << mqname << " for responding to master ..." << std::endl;
	posix_message_queue wq(mqname.c_str());
	wq.open();

	rq.notify(read_mq);

	while(true)
	{
		while(new_cmd)
		{
			//wq.send(response);
			//response.clear();
			new_cmd = false;
		}
		usleep(100000);
	}

	return 0;
}