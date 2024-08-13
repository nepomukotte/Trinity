#include <iostream>
#include <string>
#include <algorithm>
#include <algorithm>
#include <iomanip>
#include <boost/program_options.hpp>
#include <boost/asio.hpp>
#include <boost/optional.hpp>
#include <concurrency/guarded_thread.hpp>
#include <serialport/talk_to_serial.hpp>
#include "posix_message_queue.hpp"

using namespace std;
namespace ba=boost::asio;
namespace po = boost::program_options;
using fcutils::system::posix_message_queue;

std::string message = "0";
std::string response = "5050";
bool new_cmd = false;

void read_mq(sigval_t sig){
	int numRead{-1};
	// Retrive the posix message queue object from the sigval_t infos.
	posix_message_queue* mqp=static_cast<posix_message_queue*>(sig.sival_ptr);
	auto msg_size=mqp->get_msgsize();
	auto msg=new char[msg_size];
	auto mq_name=mqp->get_name();
	mqp->notify(read_mq);
	while ((numRead = mqp->receive(msg, msg_size)) > 0){
		std::cout << " Queue: " << mq_name <<", message received: " << msg << ", of size: " << numRead << std::endl;
		message = msg;
	}
	new_cmd = true;
	delete [] msg;
}

int main(int argc, char* argv[])
{
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::cout << "Staring LVPS at: " << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << std::endl;
	
	boost::optional<std::string> port;
	po::options_description desc("serial port options");
	desc.add_options()
					("help,h", "An help message")
					("port,p", po::value<boost::optional<std::string>>(&port), "The serial port to connect to");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help") || !port)
	{
		std::cout << "Pass the serial port in use like this: ./test_client_lvps -p /dev/ttyS1" << std::endl;
		return 0;
	}

	std::cout << " Using serial port: " << *port << std::endl;

	std::string mqname="/writetoPS";
	std::cout << "Creating and opening the POSIX message queue: " << mqname << " for listening to commands ..." << std::endl;
	posix_message_queue rq(mqname.c_str());
	rq.open();

	mqname="/writetoRC";
	std::cout << "Creating and opening the POSIX message queue: " << mqname << " for responding to master ..." << std::endl;
	posix_message_queue wq(mqname.c_str());
	wq.open();

	rq.notify(read_mq);

	try
	{
		// Create the listening process on the serial port
		boost::asio::io_context io_context;
		auto ttsp=fcutils::serialport::talk_to_serial<>::get(io_context,*port);
		ttsp->start_read();

		// Launch the io_context on a separate thread
		// A pointer to the member function run()
		ba::io_context::count_type (ba::io_context::*pf_run)(){&ba::io_context::run};
		// Create an launch a thread executing io_context.run()
		fcutils::concurrency::guarded_thread read_thread{pf_run,&io_context};

		std::cout << " Enter a message to write into the serial port: \n";
		while(true)
		{
			if(new_cmd)
			{
				message+='\n';
				ttsp->write(message.c_str(),message.size()+1);
				new_cmd = false;
				std::cout << " Enter a message to write into the serial port: \n";
			}
			usleep(1000);
		}
	io_context.stop();
	}
	catch (std::exception& e)
	{
		std::cerr << "Exception: " << e.what() << "\n";
	}

	return 0;
}