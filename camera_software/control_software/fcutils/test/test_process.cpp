///
///
/// @brief A new file
///
/// @file TEST_PROCESS.CPP
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 22 mar 2022 
///

#include <iostream>
#include <memory>

#include "fcutils/command/sequence.hpp"
#include "fcutils/command/sequencer.hpp"
#include "fcutils/command/sequencemanager.hpp"
#include "fcutils/command/sequencermanager.hpp"
#include "fcutils/consumer/consumer.hpp"
#include "fcutils/message/time_ht.hpp"
#include "fcutils/checker/checker.hpp"
#include "fcutils/command/sequencer.hpp"
#include "fcutils/command/sequencermanager.hpp"
#include "fcutils/command/process_common.hpp"
#include "fcutils/command/decode_process_command.hpp"
#include "fcutils/serialization/serializer.hpp"


struct executor{
	explicit executor( const std::string & s): cmd{s}{};
	executor & operator=(const std::string & s) { cmd=s; return *this;}
	void operator()(){
		std::cout <<" THis is the command: "<< cmd << std::endl;
	}
	std::string cmd;
};

template <typename Data> using time_header_t=fcutils::message::time_ht<Data, 15,fcutils::message::TIME_HT_OPTION::AT_END>;
template<typename Action> using sequencer_type=fcutils::command::sequencer<Action>;


int main(){

	using data_type = unsigned char;
	using sequence_type = fcutils::command::sequence<>;
	using manager_type = fcutils::command::seqmanager<sequence_type>;
//	using consumer_manager_type=fcutils::consumer::base_consumermanager<base_consumer,fcutils::consumer::simple_pool,time_header_t>;
	using fcutils_consumer_base_type=fcutils::consumer::consumer_base<data_type>;
	using net_consumer_manager_type=fcutils::consumer::base_consumermanager<fcutils_consumer_base_type>;
	using udp_consumer_type = fcutils::consumer::udp_consumer<data_type>;
	using tcp_consumer_type = fcutils::consumer::tcp_consumer<data_type>;
	using binaryfile_consumer_type = fcutils::consumer::binaryfile_consumer<data_type>;
	using textfile_consumer_type = fcutils::consumer::textfile_consumer<data_type,fcutils::checker::true_checker,'\n'>;
	using printbinary_consumer_type = fcutils::consumer::printbinary_consumer<data_type>;
	using printtext_consumer_type = fcutils::consumer::printtext_consumer<data_type>;
	using count_consumer_type = fcutils::consumer::count_consumer<data_type>;
	using print_stats_consumer_type = fcutils::consumer::print_stats_consumer<data_type>;
	using sequencer_base_type = fcutils::command::base_sequencer;
	using sequencermanager_type = fcutils::command::sequencermanager<sequencer_type, sequencer_base_type>;

	using process_type = fcutils::command::process_common<manager_type, sequencermanager_type, net_consumer_manager_type, std::string, data_type>;

	process_type::id_and_action_type test{static_cast<data_type>(0xFA),static_cast<data_type>(0)};
	process_type proccommon{test.id,"log"};

	fcutils::serialization::serialize_repr b;

	for( std::uint16_t i=0x50; i<0x85; ++i){
		test.action=static_cast<decltype(test.action)>(i);
		b.clear();
		b+="0x";
		b.serialize(&test);
		std::cout << " Command serialization: " << b.data() << ". \n";
		proccommon.process_msg(b.data());
	}

	fcutils::command::decode_process_command cmd_decoder;
	return 0;
}


