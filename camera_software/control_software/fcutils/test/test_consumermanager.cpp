///
///
/// @brief A new file
///
/// @file test_consumermanager.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 7 mar 2022 
///

#include <iostream>
#include <iterator>
#include <cstddef>
#include <string>

//#include "fcutils/consumer/consumermanager.hpp"
//#include "fcutils/consumer/consumer_base.hpp"
//#include "fcutils/consumer/net_consumer.hpp"
//#include "fcutils/consumer/file_consumer.hpp"
//#include "fcutils/consumer/print_consumer.hpp"
//#include "fcutils/consumer/count_consumer.hpp"
//#include "fcutils/consumer/print_stats_consumer.hpp"
//#include "fcutils/consumer/decode_consumer_command.hpp"
#include "fcutils/consumer/consumer.hpp"
#include "fcutils/message/time_ht.hpp"

struct base_consumer{
	using data_type=char;
	using data_ptr_type=char *;
	virtual void consume(const data_ptr_type, std::size_t l)=0;
	virtual void print_stats(std::ostream & o, const std::string &prefix=""){
		o << prefix << " Nothing to be printed \n";
	}
	virtual void reset(){
		std::cout << "I am resetting ... " << std::endl;
	}
};
struct print: public base_consumer{
	void consume(const  data_ptr_type d, std::size_t l){
		operator()(d,l);
	}
	void operator()(const data_ptr_type d, std::size_t l){
		std::cout << " Bytes received: " << l <<", values:\n";
		std::ostreambuf_iterator<data_type> oo (std::cout.rdbuf()); // output iterator for cout
		std::copy(d,d+l,oo);
		std::cout << std::endl;
	}
};

struct print_byte: public base_consumer{
	void consume(const data_ptr_type d, std::size_t l){
		operator()(d,l);
	}
	void operator()(const data_ptr_type d, std::size_t l){
		std::cout << " Bytes received: " << l <<", values:\n";
		for(int i=0; i<l; ++i) std::cout << std::hex<< std::to_integer<short>(static_cast<std::byte>(static_cast<unsigned char>(*(d+i)))) <<' ';
		std::cout << std::dec << std::endl;
	}
};


template <typename Data> using time_header_t=fcutils::message::time_ht<Data, 15,fcutils::message::TIME_HT_OPTION::AT_END>;
using consumer_manager_type=fcutils::consumer::base_consumermanager<base_consumer,fcutils::consumer::simple_pool,time_header_t>;
using fcutils_consumer_base_type=fcutils::consumer::consumer_base<char>;
using net_consumer_manager_type=fcutils::consumer::base_consumermanager<fcutils_consumer_base_type>;
using udp_consumer_type = fcutils::consumer::udp_consumer<char>;
using tcp_consumer_type = fcutils::consumer::tcp_consumer<char>;
using binaryfile_consumer_type = fcutils::consumer::binaryfile_consumer<char>;
using textfile_consumer_type = fcutils::consumer::textfile_consumer<char,fcutils::checker::true_checker,'\n'>;
using printbinary_consumer_type = fcutils::consumer::printbinary_consumer<char>;
using printtext_consumer_type = fcutils::consumer::printtext_consumer<char>;
using count_consumer_type = fcutils::consumer::count_consumer<char>;
using print_stats_consumer_type = fcutils::consumer::print_stats_consumer<char>;

int main(){
	std::string in;
	consumer_manager_type consumers;
	std::shared_ptr<base_consumer> ptp{new print{}};
	consumers.add_consumer("Print",ptp);
	std::shared_ptr<base_consumer> ptpb{new print_byte{}};
	consumers.add_consumer("PrintB",ptpb);
	net_consumer_manager_type net_consumers;
	std::shared_ptr<fcutils_consumer_base_type> ptn{ new fcutils::consumer::udp_consumer<char>{"127.0.0.1", 8888} };
	net_consumers.add_consumer("UDP",ptn);
	std::shared_ptr<fcutils_consumer_base_type> ptnt{ new tcp_consumer_type{"127.0.0.1",8889} };
	net_consumers.add_consumer("TCP",ptnt);
	std::shared_ptr<fcutils_consumer_base_type> ptbt{ new binaryfile_consumer_type{"Provabin",1min} };
	net_consumers.add_consumer("BIN",ptbt);
	std::shared_ptr<fcutils_consumer_base_type> pttt{ new textfile_consumer_type{"Provatext",1} };
	net_consumers.add_consumer("TXT",pttt);
	std::shared_ptr<fcutils_consumer_base_type> ptbp{ new printbinary_consumer_type{} };
	net_consumers.add_consumer("PrintB",ptbp);
	std::shared_ptr<fcutils_consumer_base_type> pttp{ new printtext_consumer_type{} };
	net_consumers.add_consumer("Print",pttp)
			("Count", std::static_pointer_cast<fcutils_consumer_base_type>(std::make_shared<count_consumer_type>(1)) )
			("ZStats", std::static_pointer_cast<fcutils_consumer_base_type>(std::make_shared<print_stats_consumer_type>(net_consumers,1)) )
			(fcutils::consumer::decode_consumer_command<char>("u 192.168.0.1:8080") )
			(fcutils::consumer::decode_consumer_command<char>("a 1") )
//			(fcutils::consumer::decode_consumer_command<char>("a 1", net_consumers))
			([](const std::string & s, auto & manager){ return fcutils::consumer::decode_consumer_command<char>(s, manager);},"a 2");
//	net_consumers.add_consumer("Count",ptcc);

	while(std::getline(std::cin,in)){
		consumers.consume_async(in.c_str(),in.size());
		net_consumers.consume_async(in.c_str(),in.size());
	}
	return 0;
}

