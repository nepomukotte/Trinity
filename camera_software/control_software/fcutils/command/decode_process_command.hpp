///
///
/// @brief Decode a message with process command
/// A process command should contain:
/// - The process id
/// - The process command
/// - Options needed by the command
/// This command coding has:
/// - Process ID: 1 Byte;
/// - Command ID: 1 Byte;
/// - Command options: Variable size, depending on the command
/// The command is supposed to be encoded into a std::string
///
/// @file decode_process_command.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 24 mar 2022 
///

#ifndef DECODE_PROCESS_COMMAND_HPP_
#define DECODE_PROCESS_COMMAND_HPP_
#include <string>
#include <iostream>
#include <type_traits>
#include <tuple>

#include "fcutils/serialization/deserializer.hpp"
#include "fcutils/utility/integer.hpp"
#include "fcutils/type_traits/type_traits.hpp"

namespace {
 namespace fctt= fcutils::type_traits;
 const std::size_t BYTE_SIZE_IN_CHAR{2};
 const std::size_t BYTE_SIZE_IN_BIT{8};
}

namespace fcutils{
namespace command{

struct decode_process_command{
using data_type = char;
using deserializer_type		= fcutils::serialization::deserialize_repr;

// The process command list. A Byte identifying the action to be executed on a process
enum class action_id: data_type {
			// Remove sequence definition
			SEQDEL=static_cast<data_type>(0x50),
			// Start sequence definition
			SEQSTART=static_cast<data_type>(0x51),
			// End sequence definition
			SEQEND=static_cast<data_type>(0x52),
			// Start sequence definition
			SEQSTARTREC=static_cast<data_type>(0x53),
			// End sequence definition
			SEQENDREC=static_cast<data_type>(0x54),
			// Execute a sequence
			SEQEXEC=static_cast<data_type>(0x55),
			// Set a low priority alarm strategy
			LALSTRATEGY=static_cast<data_type>(0x56),
			// Set an high priority alarm strategy
			HALSTRATEGY=static_cast<data_type>(0x57),
			// Reset the counters for a low priority alarm strategy
			RESETLALSTRATEGY=static_cast<data_type>(0x58),
			// Reset the counters for an high priority alarm strategy
			RESETHALSTRATEGY=static_cast<data_type>(0x59),
			// Wait
			WAIT=static_cast<data_type>(0x61),
			// Delete a low priority alarm strategy
			DELLALSTRATEGY=static_cast<data_type>(0x62),
			// Delete an high priority alarm strategy
			DELHALSTRATEGY=static_cast<data_type>(0x63),
			// Enable LOW alarm strategy
			ENABLELSTRATEGY=static_cast<data_type>(0x64),
			// Disable LOW alarm strategy
			DISABLELSTRATEGY=static_cast<data_type>(0x65),
			// Enable HIGH alarm strategy
			ENABLEHSTRATEGY=static_cast<data_type>(0x66),
			// Disable HIGH alarm strategy
			DISABLEHSTRATEGY=static_cast<data_type>(0x67),
			// Spawn a shell command
			SPAWNCMD=static_cast<data_type>(0x68),
			// Enable log file
			ENABLELOG=static_cast<data_type>(0x69),
			// Disable log file
			DISABLELOG=static_cast<data_type>(0x70),
			// Define a sequencer
			SEQUENCERREC=static_cast<data_type>(0x71),
			// Start an already registered sequencer
			SEQUENCERSTART=static_cast<data_type>(0x72),
			// Stop an already registered sequencer
			SEQUENCERSTOP=static_cast<data_type>(0x73),
			// Echo command
			ECHOCMD=static_cast<data_type>(0x74),
			// Enable log on network
			NETENABLELOG=static_cast<data_type>(0x75),
			// Disable log on network
			NETDISABLELOG=static_cast<data_type>(0x76),
			// Remove sequenceR definition
			SEQUENCERDEL=static_cast<data_type>(0x77),
			// Change SIP link
			CHANGESIPLINK=static_cast<data_type>(0x80),
			// Change SIP routing
			CHANGESIPROUTING=static_cast<data_type>(0x81),
			// Add consumer
			CONSUMERADD=static_cast<data_type>(0x82),
			// Remove consumer
			CONSUMNERDEL=static_cast<data_type>(0x83)
};

// Sizes of the option fields in Byte
enum class size_in_byte: std::size_t {
			// The process ID size
			ID = 1,
			// The action ID size
			ACTION = 1,
			// The Sequence ID size
			SEQUENCEID = 2,
			// The duration count size
			DURATION = 2,
			// The seqeuncer ID size
			SEQUENCERID = 2,
			// The time unit size
			DURATIONUNIT = 1,

};


// The ID and action
struct id_and_action{
	std::conditional_t<std::is_unsigned_v<data_type>, fcutils::uint_t<fctt::const_compare<size_in_byte,size_in_byte::ID>()*BYTE_SIZE_IN_BIT>::exact, fcutils::int_t<fctt::const_compare<size_in_byte,size_in_byte::ID>()*BYTE_SIZE_IN_BIT>::exact > id{};
	std::conditional_t<std::is_unsigned_v<data_type>, fcutils::uint_t<fctt::const_compare<size_in_byte,size_in_byte::ACTION>()*BYTE_SIZE_IN_BIT>::exact, fcutils::int_t<fctt::const_compare<size_in_byte,size_in_byte::ACTION>()*BYTE_SIZE_IN_BIT>::exact > action{};
};
using id_and_action_type 	= id_and_action;

struct sequencer_option{
	std::string id{};
	std::string seqname{};
	std::conditional_t<std::is_unsigned_v<data_type>, fcutils::uint_t<fctt::const_compare<size_in_byte,size_in_byte::DURATIONUNIT>()*BYTE_SIZE_IN_BIT>::exact, fcutils::int_t<fctt::const_compare<size_in_byte,size_in_byte::DURATIONUNIT>()*BYTE_SIZE_IN_BIT>::exact > unit{};
	std::conditional_t<std::is_unsigned_v<data_type>, fcutils::uint_t<fctt::const_compare<size_in_byte,size_in_byte::DURATION>()*BYTE_SIZE_IN_BIT>::exact, fcutils::int_t<fctt::const_compare<size_in_byte,size_in_byte::DURATION>()*BYTE_SIZE_IN_BIT>::exact > duration{};
};

using sequencer_option_type = sequencer_option;

auto decode_id_and_action(const std::string & msg){
	deserializer_type deserializer;
	// Deserialize the process ID and the actions. These are the first two bytes of the message
	id_and_action dec_id_action;
	// Deserialize the process ID and the action, use the automatic deserialization option using a representation deserializer
	auto [lastptr,val]=deserializer.deserialize(&dec_id_action,msg.c_str());
	// Check if the process_id is the same as the one in the message
	fcutils::print_vals((std::cout<<std::hex<<std::showbase),
			static_cast<std::uint32_t>(dec_id_action.id),
			static_cast<std::uint32_t>(dec_id_action.action));
	std::cout<<std::dec<<std::noshowbase;
	return lastptr;
}

auto decode_sequence_id( const std::string & msg) {
	// The sequence ID, just two bytes
	std::string seqname{msg.substr(0,fctt::const_compare<size_in_byte,size_in_byte::SEQUENCEID>()*BYTE_SIZE_IN_CHAR)};
	std::cout <<  " Sequence ID: " << seqname << std::endl;
	return seqname;
}

auto decode_duration( const std::string & msg) {
	deserializer_type deserializer;
	// an unsigned integer of the right size in bit
	fcutils::uint_t<fctt::const_compare<size_in_byte,size_in_byte::DURATION>()*BYTE_SIZE_IN_BIT>::exact duration;
	auto [lastptr,val]=deserializer.deserialize(&duration,msg.c_str());
	std::cout <<  " Duration:  " << duration << std::endl;
	return duration;
}

auto decode_sequencer( const std::string & msg) {
	deserializer_type deserializer;
	// Deserialize the sequencer ID and options.
	sequencer_option_type dec_sequencer_option;
	// Deserialize the process ID and the action, use the automatic deserialization option using a representation deserializer
	dec_sequencer_option.id=
			msg.substr(0,fctt::const_compare<size_in_byte,size_in_byte::SEQUENCERID>()*BYTE_SIZE_IN_CHAR);
	dec_sequencer_option.seqname=
			msg.substr(fctt::const_compare<size_in_byte,size_in_byte::SEQUENCERID>()*BYTE_SIZE_IN_CHAR,fctt::const_compare<size_in_byte,size_in_byte::SEQUENCEID>()*BYTE_SIZE_IN_CHAR);
	auto [lastptr,val]=deserializer.deserialize(&dec_sequencer_option.unit,msg.c_str()+(fctt::const_compare<size_in_byte,size_in_byte::SEQUENCERID>()+fctt::const_compare<size_in_byte,size_in_byte::SEQUENCEID>())*BYTE_SIZE_IN_CHAR);
	std::tie(lastptr,val) = deserializer.deserialize(&dec_sequencer_option.duration,lastptr);
	// Check if the process_id is the same as the one in the message
	fcutils::print_vals((std::cout<<std::hex<<std::showbase),
			dec_sequencer_option.id, dec_sequencer_option.seqname,
			static_cast<std::uint32_t>(dec_sequencer_option.unit),
			static_cast<std::uint32_t>(dec_sequencer_option.duration));
	std::cout<<std::dec<<std::noshowbase;
	return lastptr;
}
auto decode_sequencer_id( const std::string & msg) {
	// The sequence ID, just two bytes
	std::string sequencerid{msg.substr(0,fctt::const_compare<size_in_byte,size_in_byte::SEQUENCERID>()*BYTE_SIZE_IN_CHAR)};
	std::cout <<  " Sequence ID: " << sequencerid << std::endl;
	return sequencerid;
}

auto decode_netlog( const std::string & msg){
	return msg;
}

auto decode_consumer(const std::string & msg){
	// Consumer is decoded by the consumer manager, via the decoder object
	std::cout << " Consumer command: "<< msg << std::endl;
	return msg;
}

auto decode_consumer_id(const std::string & msg){
	// Consumer is decoded by the consumer manager, via the decoder object
	std::cout << " Consumer id: "<< msg << std::endl;
	return msg;
}

//			  if(msg.find(PROCENABLENETLOG)!=std::string::npos ) {
//				  // Duplicate the std output to an UDP address and port
//				  // extract the ip address and port from the message
//				  // Create the converter and transform the string into an array of bytes.
//				  convtochar< cmdval<proc::cmd::size::ENABLENETLOGCMD> > converter;
//				  converter.transform(msg.c_str());
//				  // Remember!!!! The position are in characters, including the leading "0x", to get it in byt, we need to divide by 2 and substract the first byte
//				  const proc::cmd::byte * address=converter.getcmd().val()+(proc::cmd::position::NETIPCMDPOS/2-1);
//				  std::string s=fcutils::network::address_to_string(address);
//				  udp_duplicators.enable(s.c_str(),std::cout);
//				  std::cout << progname << " An enable network duplicator found. Duplicating standard output on address: "<< s << std::endl;
//				    // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if(msg.find(PROCDISABLENETLOG)!=std::string::npos ) {
//				  // Disable an UDP duplicator
//				  // extract the ip address and port from the message
//				  // extract the ip address and port from the message
//				  // Create the converter and transform the string into an array of bytes.
//				  convtochar< cmdval<proc::cmd::size::ENABLENETLOGCMD> > converter;
//				  converter.transform(msg);
//				  // Remember!!!! The position are in characters, including the leading "0x", to get it in byt, we need to divide by 2 and substract the first byte
//				  const proc::cmd::byte * address=converter.getcmd().val()+(proc::cmd::position::NETIPCMDPOS/2-1);
//				  std::string s=fcutils::network::address_to_string(address);
//				  udp_duplicators.disable(s.c_str(),std::cout);
//				  std::cout << progname << " A disable network duplicator found. Disable duplication of standard output on address: "<< s << std::endl;
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if(msg.find(CONSUMERADD)!=std::string::npos ) {
//				  // Add a consumer
//				  // String is passed to a specialized function
//				  // Just strip away the header and the command ID
//				  std::string s(msg,proc::cmd::position::SPAWNCMDPOS);
//				  std::cout << progname << " A new consumer to be added, options: "<< s << std::endl;
//		    	  std::pair<std::string, consumermanager_type::consumer_ptr_type> temp=consumer::make_consumer(s);
//		    	  consumermgr.add_consumer(temp.first,temp.second);
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if(msg.find(CONSUMERDEL)!=std::string::npos ) {
//				  // Delete a consumer
//				  // Just the key is neede
//				  std::string s(msg,proc::cmd::position::SPAWNCMDPOS);
//				  std::cout << progname << " A delete consumer found for key: "<< s << std::endl;
//		    	  consumermgr.remove_consumer(s);
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }


};
} // namespace command
} // namespace fcutils




#endif // DECODE_PROCESS_COMMAND_HPP_ 
