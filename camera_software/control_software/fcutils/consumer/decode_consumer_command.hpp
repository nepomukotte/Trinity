///
///
/// @brief An helper function to decode the add consumer command string option.
/// Options string is formatted with:
/// a string specifying the label, a whitespace, a char specifying the consumer, a whitespace, the options needed to create the consumer
/// Example:
/// "u 127.0.0.1:2003", will create an UDP consumer on the ip=127.0.0.1, port=2003
/// "b /log/File 100" , will create a binary file consumer, with name prefix: /log/File, in the filename, and a logrotation of 100 minutes.
/// possible consumer and related commands are:
/// - p, for printer
/// - u, for udp
/// - t, for tcp
/// - s, for sharablememory
/// - b, for binary file
/// - x, for text file
/// - c, for counter
/// The complete list of commnands:
/// - Text printer:
/// 	w
/// - Binary printer:
/// 	p
/// - Binary file:
///		b "prefix of the filename, including root dir" "minutes for log rotation (0 for no rotation)
/// 	Example: b /home/cafagna/Blog 60 (Create a binary consumer file with logrotation and a name change after 60 minutes)
/// - Text file:
///		t "prefix of the filename, including root dir" "minutes for log rotation (0 for no rotation)
/// 	Example: x /home/cafagna/Tlog 0 (Create a text consumer file without logrotation)
/// - UDP:
///		u xxx.xxx.xxx.xxx.xxx:p
///		Where: xxx.xxx.xxx.xxx, is the IP address; p, the port number
/// 	Example: u 127.0.0.1:1234 (Create an UDP consumer on address 127.0.0.1 and port 1234)
/// - TCP:
///		t xxx.xxx.xxx.xxx.xxx:p
///		Where: xxx.xxx.xxx.xxx, is the IP address; p, the port number
/// 	Example: u 192.168.1.10:8888 (Create a TCP consumer on address 192.168.1.10 and port 8888)
/// - Sharable memory:
///		s "memory area name"
/// 	Example: s mysharedmem (Create a aharable memory consumer on memory are with name: mysharedmem)
/// - Counter:
///     c "modulo counter"
///		Example: c 100 (Create a consumer that counts the consumed packets and bytes and print the totals every 100 consumed packets)
/// - Print stats;
///     a "modulo counter"
///     Example: a 100 (Invoke the print_stats method on every consumer stored in the manager, every 100 consumed packets)
///
/// @file decode_consumer_command.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 15 mar 2022 
///

#ifndef DECODE_CONSUMER_COMMAND_HPP_
#define DECODE_CONSUMER_COMMAND_HPP_

#include <tuple>
#include <string>
#include <cctype>
#include <memory>
#include <tuple>
#include <string_view>
#include <vector>
#include <sstream>
#include <type_traits>

#include "fcutils/type_traits/type_list.hpp"
#include "fcutils/consumer/consumermanager.hpp"
#include "fcutils/consumer/consumer_base.hpp"
#include "fcutils/consumer/net_consumer.hpp"
#include "fcutils/consumer/file_consumer.hpp"
#include "fcutils/consumer/print_consumer.hpp"
#include "fcutils/consumer/count_consumer.hpp"
#include "fcutils/consumer/print_stats_consumer.hpp"
#include "fcutils/string/trim.hpp"
#include "fcutils/checker/checker.hpp"

namespace {
	struct empty_type{};
}

namespace fcutils{
namespace consumer{

enum class CONSUMER_KEY		: char 			{BINPRINTER = 'p',TXTPRINTER ='w', UDP ='u'	, TCP ='t'	, SHMEM ='s', BINFILE ='b', TXTFILE ='x', COUNTER ='c', PRINTSTATS ='a'};
enum class CONSUMER_INDEX	: std::size_t 	{BINPRINTER		 ,TXTPRINTER	 , UDP		, TCP		, SHMEM 	, BINFILE 	  , TXTFILE		, COUNTER 	  , PRINTSTATS};

template <typename Data=unsigned char, template <typename D> class Checker=fcutils::checker::true_checker>
using consumer_typelist = fcutils::type_traits::type_list<
		fcutils::consumer::printbinary_consumer<Data,Checker>,
		fcutils::consumer::printtext_consumer<Data,Checker>,
		fcutils::consumer::udp_consumer<Data,Checker>,
		fcutils::consumer::tcp_consumer<Data,Checker>,
		fcutils::consumer::printtext_consumer<Data,Checker>, // Sharable memory consumer is not yet ready, use the text printer consumer
		fcutils::consumer::binaryfile_consumer<Data,Checker>,
		fcutils::consumer::textfile_consumer<Data,Checker>,
		fcutils::consumer::count_consumer<Data,Checker>,
		fcutils::consumer::print_stats_consumer<Data,Checker>
		>;

template<typename Data=unsigned char, typename Manager=empty_type, template <typename D> class Checker= fcutils::checker::true_checker>
auto decode_consumer_command( const std::string & opt, const Manager& m=Manager{}){
	// Use a string stream .... is ugly but works....
	std::stringstream ss{opt};
	// First decode the type of consumer
	char consumer;
	ss >> consumer;
	// convert to lower case
	consumer=std::tolower(consumer);
	// STarting from C++17 we can inizialize an enum class using the undelying type ...
	CONSUMER_KEY ckey{consumer};
	// Create a typelist with the type to be used
	using consumer_types=consumer_typelist<Data,Checker>;
	// The pair to be used: Key, shared_ptr to the base type
	using consumer_base_type=fcutils::consumer::consumer_base<Data,Checker>;
	using pair_type=std::pair<std::string,std::shared_ptr< consumer_base_type > >;
	std::string stemp;
	std::string key;
	std::size_t itemp{0};
	switch(ckey){
	case CONSUMER_KEY::BINPRINTER:
	{
		constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::BINPRINTER};
		// No more options are needed
		using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
		return pair_type{"bprint", std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>())};
	}
	break;
	case CONSUMER_KEY::TXTPRINTER:
	{
		constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::TXTPRINTER};
		// No more options are needed
		using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
		return pair_type{"tprint", std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>())};
	}
	break;
	case CONSUMER_KEY::UDP:
	{
		constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::UDP};
		std::size_t colon=opt.find(':');
		if(colon==std::string::npos){
			// No port specified, temp is supposed to be a well formed ip address and the default will be used
			std::getline(ss,stemp);
		} else {
			// Get the address
			std::getline(ss,stemp,':');
			// Get the port number
			ss>>itemp;
			// Trim address, in case of leading space
			stemp=fcutils::string::trim(stemp);
		}
		std::cout <<" Creating an udp consumer with address: " << stemp <<", and port: " << itemp << std::endl;
		key="u"+stemp+"_"+std::to_string(itemp);
		using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
		return pair_type{key, std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>(stemp.c_str(),itemp))};
	}
	break;
	case CONSUMER_KEY::TCP:
	{
		constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::TCP};
		std::size_t colon=opt.find(':');
		if(colon==std::string::npos){
			// No port specified, temp is supposed to be a well formed ip address and the default will be used
			std::getline(ss,stemp);
		} else {
			// Get the address
			std::getline(ss,stemp,':');
			// Get the port number
			ss>>itemp;
			// Trim address, in case of leading space
			stemp=fcutils::string::trim(stemp);
		}
		std::cout <<" Creating a tcp consumer with address: " << stemp <<", and port: " << itemp << std::endl;
		key="t"+stemp+"_"+std::to_string(itemp);
		using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
		return pair_type{key, std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>(stemp.c_str(),itemp))};
	}
	break;
	case CONSUMER_KEY::SHMEM:
	{
		constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::SHMEM};
		// No more options are needed
		using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
		return pair_type{"shmem", std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>())};
	}
	break;
	case CONSUMER_KEY::BINFILE:
	{
		constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::BINFILE};
		ss>>stemp;
		ss>>itemp;
		std::cout <<" Creating a binary file consumer with filename: " << stemp << ", log rotate minutes: " << itemp <<std::endl;
		using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
		key="b"+stemp+"_"+std::to_string(itemp);
		return pair_type{key, std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>(stemp.c_str(),itemp))};
	}
	break;
	case CONSUMER_KEY::TXTFILE:
	{
		constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::TXTFILE};
		ss>>stemp;
		ss>>itemp;
		std::cout <<" Creating a text file consumer with filename: " << stemp << ", log rotate minutes: " << itemp <<std::endl;
		using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
		key="t"+stemp+"_"+std::to_string(itemp);
		return pair_type{key, std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>(stemp.c_str(),itemp))};
	}
	break;
	case CONSUMER_KEY::COUNTER:
	{
		constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::COUNTER};
		ss>>itemp;
		std::cout <<" Creating a count consumer with frequency: " << itemp <<std::endl;
		using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
		return pair_type{key, std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>(itemp))};
	}
	break;
	case CONSUMER_KEY::PRINTSTATS:
	{
		if constexpr(std::is_same_v<Manager,empty_type>){
			std::cout <<" ERROR!!!! No manager reference in print_stats consumer. No consumer created" <<std::endl;
		} else {
			constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::PRINTSTATS};
			ss>>itemp;
			std::cout <<" Creating a print_stats consumer with frequency: " << itemp <<std::endl;
			using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
			key="stat_"+std::to_string(itemp);
			return pair_type{key, std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>(m,itemp))};
		}
	}
	break;
	default:
		{
			std::cout <<" Ooppsss! Key: " << consumer << " not found.\n"
					<<" Using the default value: p, text printer" << std::endl;
			constexpr const CONSUMER_INDEX cidx{CONSUMER_INDEX::TXTPRINTER};
			using consumer_type=typename fcutils::type_traits::get_type<static_cast<std::size_t>(cidx),consumer_types>::type;
			return pair_type{"wrongkey", std::static_pointer_cast<consumer_base_type>(std::make_shared<consumer_type>())};
		}
		break;
	} // switch (consumer)
	// This return should have been never reached ...
	return pair_type{};
}


} // namesapce consumer
} // namespace fcutils



#endif // DECODE_CONSUMER_COMMAND_HPP_ 
