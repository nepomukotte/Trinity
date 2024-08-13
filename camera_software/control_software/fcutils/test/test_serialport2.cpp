///
///
/// @brief A new file
///
/// @file test_serialport.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 27 ott 2021 
///
#include <iostream>
#include <string>
#include <algorithm>
#include <array>
#include <variant>
#include <boost/program_options.hpp>
#include <boost/asio.hpp>
#include <boost/optional.hpp>
#include <concurrency/guarded_thread.hpp>

#include <serialport/talk_to_serial.hpp>

// Alias: boost
namespace ba=boost::asio;
namespace po = boost::program_options;
// Write your own consumer ... just copy backward into a string.

template<class Data> struct my_consumer{
	void operator()(Data* d, std::size_t l){
		// Allocate l words
		std::vector<Data> consumed(l);
		std::reverse_copy(d,d+l-1,consumed.begin());
		if(l>0){
			std::cout << " My consumer : \n"
					<< " Bytes received: " << l <<", values:\n";
			std::ostreambuf_iterator<Data> oo (std::cout.rdbuf()); // output iterator for cout
			std::copy(d,d+l,oo);
			std::cout << " Bytes consumed: " << consumed.size() <<", values:\n";
			std::copy(consumed.begin(),consumed.end(),oo);
			std::cout << "\n My consumer, job done."<< std::endl;
		}
	}
};

// Write your own checker ...
// Check if any one of a character in the array is found.
template<class Data, std::size_t Size=2> struct check_chars{
	template<typename ... Arg> check_chars(Arg... args): to_be_checked{std::forward<Arg>(args)...}{}
	size_t operator()(Data* d, size_t bytes) {
		bool found{false};
		// We read one byte at the time until we get one of the word in the to_be_checked array.
		for( auto v: to_be_checked){
			found=(*(d+bytes-1)==v)?true:found;
		}
		return found ? 0 : 1;
	}
	const std::array<Data,Size> to_be_checked;
};

using data_type=char;
template<typename Data> using print_type=fcutils::serialport::print<Data>;
template<typename Data> using my_consumer_type=my_consumer<Data>;
template<typename Data> using print_type=fcutils::serialport::print<Data>;
template<typename Data> using check_eol_type=fcutils::serialport::check_eol<Data>;
template<typename Data> using check_chars_type=check_chars<Data,2>;
using default_talk_to_serial=fcutils::serialport::talk_to_serial<>;
using my_consumer_talk_to_serial=fcutils::serialport::talk_to_serial<data_type,my_consumer_type,check_eol_type>;
using check_chars_talk_to_serial=fcutils::serialport::talk_to_serial<data_type,print_type,check_chars_type>;
using check_chars_my_consumer_talk_to_serial=fcutils::serialport::talk_to_serial<data_type,my_consumer_type,check_chars_type>;

// helper function to return the right type
template< typename MyTalkToSerial> auto get_talk_to_serial(ba::io_context &ioc, const std::string & sp){
	// Specialize the function to the types having creators, i.e. check_chars
	if constexpr (std::is_same<typename MyTalkToSerial::checker_type,check_chars_type<data_type> >::value){
		// Here I would like to init the checker with the values to be checked. So I have to pass a copy of the checker.
		return MyTalkToSerial::get(ioc,sp,check_chars_type<data_type>{'\n',';'});
	} else {
		return MyTalkToSerial::get(ioc,sp);
	}
}
int main(int ac, char* av[])
{
	// Let's look at the options: ip and port
	// A very basic usage of BOOST program options
	boost::optional<std::string> port;
	bool reverse_consumer{false},check_chars{false};
	// start adding options
	po::options_description desc("test_serialport options ");
	desc.add_options()
	        ("help,h", "An help message")
	        ("port,p", po::value<boost::optional<std::string>>(&port), "the port to connect to")
			("reverse,r", po::bool_switch(&reverse_consumer),"Use a reverse consumer (it just reverse the data stream)")
			("multiple,m", po::bool_switch(&check_chars),"Check for end of line and semicolon as end of packet markers")
	;
	// Create the option map
	po::variables_map vm;
	// Parse command line for options
	po::store(po::parse_command_line(ac, av, desc), vm);
	// Store user options into the map and the local variables
	po::notify(vm);

	if (vm.count("help")||!port) {
		std::cout << " To test the program you can create a dummy serial port using: socat \n socat -d -d pty,raw,echo=0 pty,raw,echo=0 \n"
				<< " Then run test_serialport using the ports created by socat, typically: /dev/pts/1 and /dev/pts/2. \n"
				<< " So in one terminal enter: ./test_serialport -p /dev/pts/1\n"
				<< " So in one other terminal enter: ./test_serialport -p /dev/pts/2\n\n";
		std::cout << desc << "\n";
		return 0;
	}

	std::cout << " Using serial port: " <<*port << '\n'
			<< " Using reverse consumer: " <<reverse_consumer << '\n'
			<< " Using multiple char checker: " <<check_chars << '\n';

// create a variant to hold all the sharable pointers to the talk_to_serial different types
	std::variant<default_talk_to_serial::pointer_type, my_consumer_talk_to_serial::pointer_type, check_chars_talk_to_serial::pointer_type, check_chars_my_consumer_talk_to_serial::pointer_type> ttsp;

	try
  {

    boost::asio::io_context io_context;
   	// Create the listening process on the serial port
    // fill the variant with the right type ...
    if( (!reverse_consumer)&&(!check_chars) ) ttsp=get_talk_to_serial<default_talk_to_serial>(io_context,*port);
    if( (reverse_consumer)&&(!check_chars) ) ttsp=get_talk_to_serial<my_consumer_talk_to_serial>(io_context,*port);
    if( (!reverse_consumer)&&(check_chars) ) ttsp=get_talk_to_serial<check_chars_talk_to_serial>(io_context,*port);
    if( (reverse_consumer)&&(check_chars) ) ttsp=get_talk_to_serial<check_chars_my_consumer_talk_to_serial>(io_context,*port);
    // We now need to visit the variant to execute the action
    std::visit([](auto & ptr) { ptr->start_read(); }, ttsp);

//    auto ttsp=fcutils::serialport::talk_to_serial<>::get(io_context,*port);
//    ttsp->start_read();
    // Launch the io_context on a separate thread
     // A pointer to the member function run()
 	ba::io_context::count_type (ba::io_context::*pf_run)(){&ba::io_context::run};
 	// Create an launch a thread executing io_context.run()
 	fcutils::concurrency::guarded_thread read_thread{pf_run,&io_context};
 	std::string msg;
 	std::cout << " Enter a message to write into the serial port: \n";
 	while(std::getline(std::cin,msg)){
 		// Add a carriage return
 		msg+='\n';
 		// remeber to write also the null termination ... (we need to add an extra byte)
 	    // We now need to visit the variant to execute the action
 	    std::visit([&msg](auto & ptr) { ptr->write(msg.c_str(),msg.size()+1); }, ttsp);
// 		ttsp->write(msg.c_str(),msg.size()+1);
 	 	std::cout << " Enter a message to write into the serial port: \n";
 	}
 	io_context.stop();
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception: " << e.what() << "\n";
  }

  return 0;
}



