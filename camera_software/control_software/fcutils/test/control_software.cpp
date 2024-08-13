#include <boost/asio.hpp>
#include <boost/program_options.hpp>
#include <checker/words_checker.hpp>
#include <concurrency/sync_queue.hpp>
#include <networking/client.hpp>
#include <networking/server.hpp>
#include <networking/session.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <string_view>
#include <chrono>
#include <ctime>
#include <cstring>
#include <thread>
#include <cstddef>

#include "string/char_to_hex.hpp"
#include "system/posix_message_queue.hpp"

#define POWER_ON_SEQ	0x01
#define INIT_SEQ		0x02
#define CONFIGURE_SEQ	0x03
#define START_SEQ		0x04
#define STOP_SEQ		0x05
#define POWER_OFF_SEQ	0x06

std::string SEQ_DIR = "/home/cherenkov/Programs/control_software/sequences/";

// Alias: std
namespace ch = std::chrono;
// Alias: boost
namespace ba=boost::asio;
namespace po = boost::program_options;

namespace fcc=fcutils::checker;
namespace fcn=fcutils::networking;
namespace fcs=fcutils::string;
namespace fcsy=fcutils::system;

// A message queue
using my_queue=fcutils::concurrency::sync_queue<std::string,false>;
using my_posix_queue=fcutils::system::posix_message_queue;

std::string mqname[2] = {"/CStoRC","/RCtoCS"};
my_posix_queue wq(mqname[0].c_str());	// the posix message to write to
my_posix_queue rq(mqname[1].c_str());	// the posix message to read from

const std::size_t GCCHEADER_SIZE{8},GCCHEADER_CHAR_SIZE{GCCHEADER_SIZE*2},GCCSTATUS_MAX_SIZE{104};

struct message_consumer{
	my_queue & head, & user;
	int i{0};
	const char eotx{0x04};

	message_consumer(my_queue & q, my_queue & u) :head{q}, user{u} {};
	// To read the TCP client or server sessions
	void operator()(){
		++i;
		std::string msg;
		head.get(msg);
		if(msg[0]==eotx) {
			std::cout << " Consumer end of transmission detected ... \n";
			user.put(std::move(msg));
			return;
		}
		// Strip away the GCC header
		std::string_view sv{msg};
		if( msg.size()>GCCHEADER_CHAR_SIZE ) sv.remove_prefix(GCCHEADER_CHAR_SIZE);
		std::cout << "Consumer i: "<< i <<", full GCC message: 0x"<< msg <<", payload: 0x"<< sv<< std::endl;
		// Transmit the message to the user
		user.put(std::move(msg));
		// post back this member function into the thread and wait for another message
		(*this)();
	};

};

template <typename Data> struct message_producer{
	using data_type=Data;
	my_queue & tail;
	message_producer(my_queue & q) :tail(q) {};
	void operator()(data_type * d, std::size_t length){
		std::string msg{fcutils::string::char_to_hex(d,length,"")};
		std::cout << " A message received: " << msg <<'\n';
		tail.put(std::move(msg));
		return;
	};
	const char eotx{0x04};
};

// We need here the message queues to guarantee usage of a standard function to the posix message queue ...
// To be solved somehow in the future ...

my_queue uq;

// To read the posix message queue
void read_posix_mq(sigval_t sig){
	int numRead{-1};
	// Retrive the posix message queue object from the sigval_t infos.
	my_posix_queue* mqp=static_cast<my_posix_queue*>(sig.sival_ptr);
	auto msg_size=mqp->get_msgsize();
	auto msg=new char[msg_size];
	auto mq_name=mqp->get_name();
	mqp->notify(read_posix_mq);
	numRead = mqp->receive(msg, msg_size);
	if(numRead > 0){
		std::cout << " Queue: " << mq_name <<", message received: " << msg << ", of size: " << numRead << std::endl;
		std::string temp;
		temp.reserve(numRead);
		//temp="CT";
		temp.append(msg,numRead);
		uq.put(std::move(temp));
	}
	delete [] msg;
}

using data_type=char;
using key_type=char;
using checker_return_type=std::size_t;
const bool DEBUG{false}; // Set this to true for more DEBUG printout from the checker
template<typename Data> using gcc_consumer_type=message_producer<Data>;
template<typename Data> using gcc_checker_type=fcc::words_checker<GCCHEADER_SIZE,checker_return_type,fcc::decode_header,Data,key_type,DEBUG>;
template<typename Protocol, bool Echo> using gcc_session_type=fcn::session<Protocol, Echo, gcc_consumer_type, gcc_checker_type >;
// The client using the checker for a GCC header
using gcc_client_type=fcn::base_client< fcn::tcp, false, gcc_session_type >;

std::string find_seq_name(const std::string& seq_msg)
{
	std::string seq_str = seq_msg.substr(0, 2);
	int seq_id = stoi(seq_str,0,16);
	std::cout << "Sequence ID: " << seq_id << " was received." << std::endl;
	std::string filename = "";
	switch(seq_id)
	{
		case POWER_ON_SEQ:
			filename = SEQ_DIR + "power_on_seq.txt";
			break;
		case INIT_SEQ:
			filename = SEQ_DIR + "init_seq.txt";
			break;
		case CONFIGURE_SEQ:
			filename = SEQ_DIR + "config_seq.txt";
			break;
		case START_SEQ:
			filename = SEQ_DIR+ "start_daq_seq.txt";
			break;
		case STOP_SEQ:
			filename = SEQ_DIR + "stop_daq_seq.txt";
			break;
		case POWER_OFF_SEQ:
			filename = SEQ_DIR + "power_off_seq.txt";
			break;
		default:
			break;
	}
	return filename;
}

int main(int ac, char* av[])
{

	// Let's look at the options: ip and port
	// A very basic usage of BOOST program options
	std::string address,write_mq,read_mq;
	ba::ip::port_type port;
	int delay{0};
	bool timestamp{false},nostatus{false};
	// start adding options
    po::options_description desc("control software options ");
    desc.add_options()
        ("help,h", "An help message")
        ("ip", po::value<std::string>(&address)->default_value("192.168.2.1"), "the IP address to talk to")
        ("port,p", po::value<ba::ip::port_type>(&port)->default_value(20020U), "the port to connect to")
		("delay,d", po::value<int>(&delay)->default_value(30), "the delay between broadcast (seconds)")
		("timestamp,t", po::bool_switch(&timestamp),"use timestamp instead of a counter")
		("nostatus,n", po::bool_switch(&nostatus),"turn off automatic status messages")
		//("read_queue,r", po::value<std::string>(&read_mq)->default_value("/RCtoCS"), "the posix message to read from")
		//("write_queue,w", po::value<std::string>(&write_mq)->default_value("/CStoRC"), "the posix message to write into")
    ;
	// Create the option map
    po::variables_map vm;
    // Parse command line for options
    po::store(po::parse_command_line(ac, av, desc), vm);
    // Store user options into the map and the local variables
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    std::cout << " IP address: " << address <<":"<<port<<", delay (sec): "<<delay<< '\n';
    std::string s=(timestamp)?"timestamp":"message counter";
    std::cout << "Using "<<s<<" as status message \n";
    std::cout << "Reading from a posix message_queue: "<< mqname[1].c_str() <<'\n';
    std::cout << "Writing into a posix message_queue: "<< mqname[0].c_str() <<'\n';

    // Here the checker for words of size char, storable in a checker map with a char key.
	// Here the number of words to be checked is the size of the header
    // This object will be copied into a client or a server.
	gcc_checker_type<data_type> gcc_checker{
		'g',
		{ // The header byte constant values to be checked for only 7 of 8 bytes will be used,
			static_cast<data_type>(0xAA), // The header ID, four bytes 0xAA5555AA
			static_cast<data_type>(0x55),
			static_cast<data_type>(0x55),
			static_cast<data_type>(0xAA),
			static_cast<data_type>(0x00), // The message type: 0x0, for received msg; 0x1 for status, i.e. transmitted, message
			static_cast<data_type>(0x00), // Mandatory empty bytes, two bytes
			static_cast<data_type>(0x00)
		},
		fcc::decode_header<bool,std::size_t>{8,8} // The number of words to be checked upon extraction of the size from the offset.
	}; // end of words_checker constructor

	// Create the message queues
	//my_posix_queue rq(read_mq.c_str());
	rq.open();
	//my_posix_queue wq(write_mq.c_str());
	wq.open();

    // Try to deal with the TCP stream objects
    try{
    	// start with an ASIO context
    	ba::io_context ioc;
    	// We need a consumer.
    	// The consumer is supposed to be passed by reference.
    	// This is because a consumer can be also a container that could be filled with new consumers at run time.
    	// For this example we use a syncronized message queue producer that will serialize the message into a string, and put into the queue.
    	my_queue mq;
		gcc_client_type::session_type::consumer_type sq_producer{mq};
    	// Then we need a message receiver
    	message_consumer sq_consumer{mq,uq};
    	// Start the thread with the message queue receiver
    	rq.notify(read_posix_mq);
    	// Start the message receiver in a separate thread
    	std::thread my_receiver{sq_consumer};
    	// Construct a signal set registered for process termination.
    	ba::signal_set signals(ioc, SIGINT, SIGTERM);
    	// Start an asynchronous wait for one of the signals to occur.
    	signals.async_wait([&mq](const boost::system::error_code& ec, int signal_number){
    		std::cout << " Signal received: " << signal_number <<'\n';
    		if(!ec) {
    			std::cout << " Signal handler. No error \n";
    		} else {
    			std::cout << " signal handler. Error : " << ec.message() << std::endl;
    		}
    		// Anyway, if we receive a signal means that something went wrong. So let's stop the message queue
    		const char eotx{0x04};
    		mq.put(my_queue::message_type{eotx});
    	});
    	// The client TCP to be used with a session checking for the GCC header
    	gcc_client_type gcc_client{ioc,sq_producer,gcc_checker};
    	auto ep=fcn::tcp::endpoint{boost::asio::ip::make_address(address),port};
    	// Add the address and port and try to connect
		// gcc_client.async_connect(boost::asio::ip::make_address(address),port);
    	// Set the reuseaddress and reuseports options to the socket
    	gcc_client.async_connect(ep);
    	gcc_client.set_option(ep,boost::asio::socket_base::reuse_address{true});
    	gcc_client.set_option(ep,boost::asio::detail::socket_option::boolean<SOL_SOCKET, SO_REUSEPORT>{true});
    	// Start the client
    	gcc_client.start();
    	// Start the context on another thread
        // A pointer to the member function run(). Is needed because is an oveloaded function
        ba::io_context::count_type (ba::io_context::*pf_run)(){&ba::io_context::run};
        // Create an launch a thread executing io_context.run()
        std::thread t{pf_run,&ioc};
       	// create a time point converted to time_t (C style) for an easier text conversion
       	std::time_t now=ch::system_clock::to_time_t(ch::system_clock::now());
       	// create a time duration for the wait state
       	auto wd=ch::seconds{delay};
       	int count{0};
       	// Prepare a standard header needed by the GCC message
       	const std::size_t GCCHEADER_SIZE{8},GCCSTATUS_MAX_SIZE{104};
       	unsigned char gccheader[GCCHEADER_SIZE]={static_cast<unsigned char>(0xAA),
       			static_cast<unsigned char>(0x55),static_cast<unsigned char>(0x55),
				static_cast<unsigned char>(0xAA),static_cast<unsigned char>(0x01),
				static_cast<unsigned char>(0x00),static_cast<unsigned char>(0x00) };
       	// Just the null termination char
       	char gcctrailer='\0';
       	// create the buffer sequence
       	std::array<ba::const_buffer,2> status_msg;
       	// A pointer to the message
       	char *cnow{nullptr};
       	// The message length
       	std::size_t length{0};
       	std::string msg;
       	// End of Transmission according to ASCII coding
		const char eotx{0x04};
       	while(1){
       		if(!nostatus){
				++count;
				std::cout<< "Sending status message no. "<<count<<" ...\n";
				if(timestamp){
					// this is in principle a static array of char ...
					auto now=ch::system_clock::to_time_t(ch::system_clock::now());
					cnow=std::ctime(&now);
					// Remember, std::strlen returns the lenght without the null termination
					length=cnow?std::strlen(cnow)+1:0;
					// Set maximum length to 146 bytes
					length=(length<=GCCSTATUS_MAX_SIZE)?length:GCCSTATUS_MAX_SIZE;
				} else {
					cnow=reinterpret_cast<char *>(&count);
					length=sizeof(count);
				}
				// Store the message length into the gcc header
				gccheader[7]=static_cast<char>(length);
	//       		stream<<std::ctime(&now);
	//       		stream << count;
				std::cout << " Writing the message :"<< cnow << ", of length: "<< length <<std::endl;
				// Create the buffer sequence
				status_msg[0]=ba::buffer(gccheader,GCCHEADER_SIZE);
				status_msg[1]=ba::buffer(cnow,length);
				gcc_client.write(status_msg);
				// Put this thread on wait for 30 seconds
	//           	std::this_thread::sleep_for(wd);
				// Instead of sleeping, just check if there is any message in the queue
       		} // if(!nostatus)
       		msg.clear();
       		uq.get(msg,wd);
    		if(msg[0]==eotx) {
    			std::cout << " End of transmission detected ... That's all folks!!!\n";
    			return 0;
    		}
    		// Here we get a timeout or a message was set.
    		if(!msg.empty()) {
    			if( msg[0]=='C' && msg[1]=='T'){
    				// Here we got a message from the CT forward it to the GCC
    				// Store the message length into the gcc header
    				length=msg.size();
    				gccheader[7]=static_cast<char>(length);
    				std::cout << " Writing the message :"<< msg << ", of length: "<< length <<std::endl;
    				// Create the buffer sequence
    				status_msg[0]=ba::buffer(gccheader,GCCHEADER_SIZE);
    				status_msg[1]=ba::buffer(msg.c_str(),length);
    				gcc_client.write(status_msg);
    			} else {
    				// Here we got a message from the GCC
    				// Strip away the GCC header
    				std::string_view sv{msg};
    				if(msg.size()>GCCHEADER_CHAR_SIZE) sv.remove_prefix(GCCHEADER_CHAR_SIZE);
    				std::cout << " A message from GCC. Full GCC message: 0x"<<msg <<", payload: 0x"<< sv<< std::endl;
    				// Transmit the message to the user
    				std::string umsg{sv};
    				// Checking to see if the message is a sequence or a single command
    				std::string umsg_id_str = umsg.substr(0, 2);
					int umsg_id = stoi(umsg_id_str,0,16);
					if(umsg_id == 0xAA)
					{
						std::string seq_payload = umsg.substr(2, 16);
        				std::cout << "A sequence was received. Starting to execute the sequence. Here is the sequence payload: " << seq_payload << '\n';
						std::string filename = find_seq_name(seq_payload);

						std::string line;
						std::cout << "Loading the filename: " << filename << std::endl;
						std::ifstream infile;
						infile.open(filename.c_str());
						while(std::getline(infile, line))
						{
							std::getline(infile, line);
							std::istringstream iss(line);
							std::string command = iss.str();
							std::cout << "sequence message: " << command << std::endl;
							if(command.substr(0,1) != "*"){
								wq.send(command.c_str(), command.size());								
							}
							else{
								std::cout << "skipping this command ..." << std::endl;
							}
							int delay = 0;
							std::getline(infile, line);
							std::istringstream iss_delay(line);
							iss_delay >> delay;
							std::cout << "delaying for " << delay << " seconds."<< std::endl;
							sleep(delay);
						}
						std::cout << "Finished loading this sequence: " << filename << std::endl;
						infile.close();
					} else{
						std::cout << "A message was received and sent to the message queue: " << umsg << '\n';
    					wq.send(umsg.c_str(),umsg.size());
					}
    			}
    		} // if(!msg.empty)
       	} // while(1)
    } // try
    // Catch ASIO exceptions here
    catch (std::exception &e)
        {
            std::cerr << "Exception: " << e.what() << "\n";
        }//

	return 0;
}

