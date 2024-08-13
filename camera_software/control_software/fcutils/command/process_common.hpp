///
///
/// @brief Common process utilities. A struct that handles a sequence and sequencer managers along with logger utilities.
///
/// @file process_common.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 21 mar 2022 
///

#ifndef PROCESS_COMMON_HPP_
#define PROCESS_COMMON_HPP_
#include <string>
#include <chrono>
#include <iosfwd>
#include <type_traits>
#include <thread>
#include <errno.h> // for program_invocation_short_name
#include <string_view>

#include "fcutils/chrono/chrono.hpp"
#include "fcutils/string/program_name_and_date.hpp"
#include "fcutils/serialization/deserializer.hpp"
#include "fcutils/utility/print_vals.hpp"
#include "fcutils/iostream/iostream.hpp"

namespace {
	using namespace std::chrono_literals;
	namespace proc{
	namespace cmd {
		// Command sizes
		typedef unsigned char byte;
		// The trailer
		const unsigned char EOC=0xEC;
		// All sequencer, sequences, anything is 2 bytes long
		const std::string ALL("0000");
		namespace size{
			// The PROC Command size, in characters  including the leading 0x characters
			const std::size_t PROCCMDCHARSIZE=6;
			// The PROC Header size, in characters  including the leading 0x characters
			const std::size_t PROCHEADCHARSIZE=4;
			// The sequence name size, in bytes
			const std::size_t SEQIDSIZE=2;
			// The sequencer name size, in bytes
			const std::size_t SEQUENCERIDSIZE=2;
			// The time or counter unit size, in bytes
			const std::size_t TIMECNTUNITSIZE=1;
			// The time or counter value size, in bytes
			const std::size_t TIMECNTVALSIZE=2;
			// The Alarm mask size, in bytes
			const std::size_t LALARMMASKSIZE=5;
			const std::size_t HALARMMASKSIZE=1;
			// The Alarm action counter size, in bytes
			const std::size_t ALARMCOUNTERSIZE=1;
			// The header, command, CRC and trailer size, in bytes
			const std::size_t CMDADDRESSSIZE=1;
			const std::size_t CMDACTIONSIZE=1;
			const std::size_t CMDCHANNELSIZE=1;
			const std::size_t CMDCRCSIZE=2;
			const std::size_t CMDEOCSIZE=1;
			const std::size_t CMDDATASIZE=4;
			const std::size_t DATA=2;
			// The SIP link and routing sizes
			const std::size_t SIPLINKSIZE=1;
			const std::size_t SIPROUTINGSIZE=1;
			// The IP address and port size for UDP log.
			// Four (4) bytes are for the IP address 2 for the port
			const std::size_t IPADDRESSSIZE=7;
			const std::size_t ASKCMD=CMDADDRESSSIZE+CMDACTIONSIZE+CMDCHANNELSIZE+CMDEOCSIZE+CMDCRCSIZE; // It should be 6
			const std::size_t REPLYCMD=8;// It should be 8
			const std::size_t SETCMD=7;// It should be 7
			const std::size_t OPCMD=5;// It should be 5
			const std::size_t LALARMSTRATEGYCMD=CMDADDRESSSIZE+CMDACTIONSIZE
					+LALARMMASKSIZE+SEQIDSIZE+ALARMCOUNTERSIZE+CMDEOCSIZE+CMDCRCSIZE; // It should be 14
			const std::size_t HALARMSTRATEGYCMD=CMDADDRESSSIZE+CMDACTIONSIZE
					+HALARMMASKSIZE+SEQIDSIZE+ALARMCOUNTERSIZE+CMDEOCSIZE+CMDCRCSIZE; // It should be 10
			const std::size_t LALARMRESETSTRATEGYCMD=CMDADDRESSSIZE+CMDACTIONSIZE
					+LALARMMASKSIZE+CMDEOCSIZE+CMDCRCSIZE; // It should be 10
			const std::size_t HALARMRESETSTRATEGYCMD=CMDADDRESSSIZE+CMDACTIONSIZE
					+HALARMMASKSIZE+CMDEOCSIZE+CMDCRCSIZE; // It should be 6
			const std::size_t SEQUENCERCMD=CMDADDRESSSIZE+CMDACTIONSIZE+SEQUENCERIDSIZE
					+SEQIDSIZE+TIMECNTUNITSIZE+TIMECNTVALSIZE+CMDEOCSIZE+CMDCRCSIZE; // It should be 12
			const std::size_t ENABLENETLOGCMD=CMDADDRESSSIZE+CMDACTIONSIZE+IPADDRESSSIZE+CMDEOCSIZE+CMDCRCSIZE; // it should be 12
			const std::size_t SETSIPLINKCMD=CMDADDRESSSIZE+CMDACTIONSIZE+SIPLINKSIZE+CMDEOCSIZE+CMDCRCSIZE; // it should be 12
			const std::size_t SETSIPROUTECMD=CMDADDRESSSIZE+CMDACTIONSIZE+SIPROUTINGSIZE+CMDEOCSIZE+CMDCRCSIZE; // it should be 12
		} // namespace size
	}// namesapce cmd
	} // namespace proc
}

namespace fcutils{
namespace command{

template<typename SequenceManager, typename SequencerManager, typename ConsumerManager, typename Logger=fcutils::iostream::ostream_logger, typename Data= unsigned char,
		typename Duration=std::chrono::milliseconds>
struct process_common {
	using self_type						= process_common; //<SequenceManager, SequencerManager, Logger, Data, Duration>;
	using data_type						= Data;
	using sequence_manager_type			= SequenceManager;
	using sequence_id_type				= typename sequence_manager_type::id_type;
	using sequencer_manager_type		= SequencerManager;
	using sequencer_id_type				= typename sequencer_manager_type::id_type;
	using sequencer_ptr_type			= typename sequencer_manager_type::sequencer_ptr_type;
	using consumer_manager_type			= ConsumerManager;
	using logger_type					= Logger;
	using duration_type					= Duration;
	using deserializer_type				= fcutils::serialization::deserialize_repr;
	using net_duplicator_manager_type	= fcutils::iostream::net_duplicator_manager<>;
	struct id_and_action{
		data_type id{};
		data_type action{};
	};
	using id_and_action_type 	= id_and_action;

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

	process_common(data_type procid, const std::string & loggername, const std::string & progname=program_invocation_short_name):
			procid_(procid), proc_byte_helper(procid), progname(progname), logger(loggername), seqdef_(false)
		{
//			PROCHEADER=proc_byte_helper_.str(procid);
//			SEQDEL=proc_byte_helper_.str(proc::cmd::action::SEQDEL);
//			SEQSTART=proc_byte_helper_.str(proc::cmd::action::SEQSTART) ;
//			SEQEND=proc_byte_helper_.str(proc::cmd::action::SEQEND) ;
//			SEQSTARTREC=proc_byte_helper_.str(proc::cmd::action::SEQSTARTREC) ;
//			SEQENDREC=proc_byte_helper_.str(proc::cmd::action::SEQENDREC);
//			SEQEXEC=proc_byte_helper_.str(proc::cmd::action::SEQEXEC);
//			WAIT=proc_byte_helper_.str(proc::cmd::action::WAIT);
//			ENABLELOG=proc_byte_helper_.str(proc::cmd::action::ENABLELOG);
//			DISABLELOG=proc_byte_helper_.str(proc::cmd::action::DISABLELOG);
//			OPCMDHEADER=lower_hkop.str();
//			SSCMDHEADER=lower_hksys.str();
//			SEQUENCERREC=proc_byte_helper_.str(proc::cmd::action::SEQUENCERREC);
//			SEQUENCERSTART=proc_byte_helper_.str(proc::cmd::action::SEQUENCERSTART);
//			SEQUENCERSTOP=proc_byte_helper_.str(proc::cmd::action::SEQUENCERSTOP);
//			SEQUENCERDEL=proc_byte_helper_.str(proc::cmd::action::SEQUENCERDEL);
//			PROCSPAWNCMD=proc_byte_helper_.str(proc::cmd::action::SPAWNCMD);
//			PROCECHOCMD=proc_byte_helper_.str(proc::cmd::action::ECHOCMD);
//			PROCENABLENETLOG=proc_byte_helper_.str(proc::cmd::action::NETENABLELOG);
//			PROCDISABLENETLOG=proc_byte_helper_.str(proc::cmd::action::NETDISABLELOG);
//			CONSUMERADD=proc_byte_helper_.str(proc::cmd::action::CONSUMERADD);
//			CONSUMERDEL=proc_byte_helper_.str(proc::cmd::action::CONSUMNERDEL);
		}


// Data Member
	data_type proc_byte_helper;
	fcutils::string::program_name_and_date progname;
	net_duplicator_manager_type net_duplicators;
	// create the logger and immediatelly enable log in case is needed
	logger_type logger;
	sequence_manager_type sequencemgr;
	//// Define the sequencer manager type
	sequencer_manager_type sequencermgr;
	/// Add a consumer manager
	consumer_manager_type consumermgr;
private:
	// The private data part
	data_type procid_;
	// Sequence recording flag. REMEMBER it must be initialized always to false
	bool seqdef_{false};
	// Last Sequence name to be recorded
	std::string seqname_;
	// A flag to enable command execution. Is true if the command is not a command for the process
	bool cmdexec_{false};
	// A deserializer to extract data field from the command string
	deserializer_type deserializer_;
	// tHE NAME OF THE ANONYMOUS SEQUENCE
	const std::string ANOMSEQ_{"noname"};

public:
	// Member functions

	template<typename Ret> Ret get_id() const { return static_cast<Ret>(procid_); }

	template <typename D> void sleep_for(const D & d){
		// If Duration is an integer then create a standard duration, otherwise cast to the standard duration
		if constexpr (std::is_integral<D>::value){
			duration_type wait{d};
			std::cout << progname << "Waiting for "<< wait.count() << fcutils::chrono::suffix<duration_type>() << std::endl;
			std::this_thread::sleep_for(wait);
		} else {
			duration_type wait{std::chrono::duration_cast<duration_type>(d)};
			std::cout << progname << "Waiting for "<< wait.count() << fcutils::chrono::suffix<duration_type>() << std::endl;
			std::this_thread::sleep_for(wait);
		}
	}

	// we need a copy here ... the function should be used for callback ....
	void spawn_cmd(const std::string &s){
		//Launch child process
		std::cout << progname << " Spawing command: " << s << std::endl;
		if(0 != std::system(s.c_str()))
			std::cout << progname << " ERROR ->>> Error spawing command: " << s << std::endl;
	}


	void delete_sequence( const std::string & seqname){
			// a sequence DELETE is found.
			seqdef_=false;
//			// Extract the sequence name (ID) from the message and remove it
//			// Remember the string first two char are the hex identifier "0x", then there is the address, the command and the ID
//			std::size_t start=2+2+2;
//			auto seqname=s.substr(0,proc::cmd::size::SEQIDSIZE*2 );
			// reset the command size to avoid execution to the serial
			auto cmdsize=0;
			std::cout << progname << " A delete sequence found. Sequence ID: " << seqname << std::endl;
			// remove the sequence
			sequencemgr.remove_sequence(seqname);
	}
	//			if( ( dec_id_action!=std::string::npos )  ){
	//				// a sequence DELETE is found.
	//				seqdef_=false;
	//				// Extract the sequence name (ID) from the message and remove it
	//				// Remember the string first two char are the hex identifier "0x", then there is the address, the command and the ID
	//				std::size_t start=2+2+2;
	//				seqname=msg.substr(start,proc::cmd::size::SEQIDSIZE*2 );
	//				// reset the command size to avoid execution to the serial
	//				cmdsize=0;
	//				std::cout << progname << " A delete sequence found. Sequence ID: " << seqname << std::endl;
	//				// remove the sequence
	//				sequencemgr.remove_sequence(seqname);
	//			}

	void end_sequence ( const std::string & seqname){
		// a sequence definition is over, just connect to the signal
		  seqdef_=false;
		  std::cout << progname << " An end of sequence found. Executing the sequence ..." << std::endl;
		  // execute the sequence
		  sequencemgr.execute_sequence(ANOMSEQ_);
		  sequencemgr.remove_sequence(ANOMSEQ_);
		  // reset the command size to avoid execution to the serial
		  auto cmdsize=0;
	}

	void startrec_sequence ( const typename sequence_manager_type::id_type & seqname){
		  // Start a sequence recording
		  seqdef_=true;
		  // reset the command size to avoid execution to the serial
		  auto cmdsize=0;
		  std::cout << progname << " A start of sequence recording found. Sequence ID: " << seqname << std::endl;
		  // Add a new empty sequence
		  add_sequence(seqname);
	}

	void endrec_sequence () {
	// an end of recording sequence found. Just flag the end of a sequence definition
			seqdef_=false;
			std::cout << progname << " An end of sequence recording found. " << std::endl;
			// reset the command size to avoid execution to the serial
		auto cmdsize=0;
	}

	void add_sequence(const typename sequence_manager_type::id_type & id ){
		  std::cout << progname << " A create a new sequence found. Sequence ID: " << id << std::endl;
		  // Add a new empty sequence
		  sequencemgr.add_sequence(id);

	}

	void add_cmd_to_sequence(const typename sequence_manager_type::id_type & id, const typename sequence_manager_type::command_type & cmd) {
		  std::cout << progname << " An add new command found.  Sequence ID: " << id << std::endl;
		  // Add a new empty sequence
		  sequencemgr.add_cmd_to_sequence(id,cmd);

	}
	void exec_sequence( const std::string & seqname){
		// execute a sequence
		sequencemgr.execute_sequence(seqname);
		// reset the command size to avoid execution to the serial
		auto cmdsize=0;

	}

	void enable_log( std::ostream & o = std::cout) {
		  // Enable a log file of the std output
		// Open a log file with the ID and time stamps
		logger.enable(o);
		std::cout << progname << " An enable log found. Duplicating standard output on log file: "<< logger.get_name() << std::endl;
		// reset the command size to avoid execution to the serial
		auto cmdsize=0;
	}

	void disable_log( std::ostream & o = std::cout) {
		// Disable a log file of the std output
		logger.disable(o);
		std::cout << progname << " A disable log found. Disable duplication of standard output on log file: "<< logger.get_name() << std::endl;
		// reset the command size to avoid execution to the serial
		auto cmdsize=0;
	}

	template<class D, class Action > void add_sequencer(const sequencer_id_type & id, const D &val, std::shared_ptr<Action> a  ){
		  std::cout << progname << " A record a new sequencer found. ID (hex): "
				  << id
				  << ", Time value: " << val <<  std::endl;
		  sequencermgr.add_sequencer(id,val,a);
		  // reset the command size to avoid execution to the serial
		  auto cmdsize=0;
	}

	template<class D>
	void add_sequencer(const sequencer_id_type & id, const D &val,
			const sequence_id_type & seqname  ){
		  std::cout << progname << " A record a new sequencer found. ID (hex): "
				  << id
				  << ", Time value: " << val <<  std::endl;
		  sequencermgr.add_sequencer(id,val,sequencemgr[seqname]);
		  // reset the command size to avoid execution to the serial
		  auto cmdsize=0;
	}

	void start_sequencer(const sequencer_id_type & id){
		  std::cout << progname << " A start sequencer found. Sequencer ID (hex): " << id << std::endl;
		  sequencermgr.start_sequencer(id);
		  auto cmdsize=0;
	}

	void start_all_sequencer(){
		  std::cout << progname << " A start all sequencer found" << std::endl;
		  sequencermgr.start_all_sequencer();
		  auto cmdsize=0;
	}

	void stop_sequencer(const sequencer_id_type & id){
		  std::cout << progname << " A stop sequencer found. Sequencer ID (hex): " << id << std::endl;
		  sequencermgr.stop_sequencer(id);
		  auto cmdsize=0;
	}

	void stop_all_sequencer(){
		  std::cout << progname << " A stop all sequencer found" << std::endl;
		  sequencermgr.stop_all_sequencer();
		  auto cmdsize=0;
	}

	void delete_sequencer(const sequencer_id_type & id){
		  std::cout << progname << " A delete sequencer found. Sequencer ID (hex): " << id << std::endl;
		  // First stop the sequencer
		  sequencermgr.stop_sequencer(id);
		  // Then remove it
		  sequencermgr.remove_sequencer(id);
		  auto cmdsize=0;
	}

	void delete_all_sequencer(){
		  std::cout << progname << " A delete all sequencer found" << std::endl;
		  // First stop all sequencer
		  sequencermgr.stop_all_sequencer();
		  // Then remove them all
		  sequencermgr.remoce_all_sequencer();
		  auto cmdsize=0;
	}

	std::string process_echo( std::string & s){
		// Just echo the msg received, as it is, adding the program name
		std::cout << progname << " Echo command. " << std::endl;
		// reset the command size to avoid dispatch to the serial port
		auto cmdsize=0;
		std::string temp{s};
		temp+=progname.get_name();
		return temp;

	}

	void enable_netlog(const std::string & cmd, std::ostream & o=std::cout){
		// Duplicate the std output to a net protocol, address and port
		// Is the net_duplicator_manager responsibility to decode the command, that is the transport protocol: UDP or TCP, the IP address and the port.
		net_duplicators.enable(cmd,o);
		std::cout << progname << " An enable network duplicator found. Duplicating standard output on net duplicator: "<< cmd << std::endl;
		// reset the command size to avoid execution to the serial
		auto cmdsize=0;
	}
	void disable_netlog( const std::string & cmd, std::ostream & o = std::cout){
		// Disable a net duplicator.
		// Is the net_duplicator_manager responsibility to decode the command, that is the transport protocol: UDP or TCP, the IP address and the port.
		net_duplicators.disable(cmd,o);
		std::cout << progname << " A disable network duplicator found. Disable net duplication of standard output: "<< cmd << std::endl;
		// reset the command size to avoid execution to the serial
		auto cmdsize=0;
	}

	template <typename Decoder>
	void add_consumer(Decoder && decoder, const std::string & cmd){
		// Add a consumer
		// Consumer is decoded by the consumer manager, via the decoder object
		std::cout << progname << " A new consumer to be added: "<< cmd << std::endl;
		consumermgr.add_consumer(std::forward<Decoder>(decoder),cmd);
		// reset the command size to avoid execution to the serial
		auto cmdsize=0;
	  }

	template<typename Ptr=typename consumer_manager_type::consumer_ptr_type>
	void add_consumer(const typename consumer_manager_type::key_type & id, Ptr && c){
		// Add a consumer pointed by c with id
		std::cout << progname << " A new consumer to be added with id: "<< id << std::endl;
		consumermgr.add_consumer(id,std::forward<Ptr>(c));
		// reset the command size to avoid execution to the serial
		auto cmdsize=0;
	}

	template<typename Node=typename consumer_manager_type::consumer_map_type::node_type>
	void add_consumer(Node && n){
		// Add a consumer pointed by c with id
		std::cout << progname << " A new node in the consumer map to be added with id: "<< n.first << std::endl;
		consumermgr.add_consumer(std::forward<Node>(n));
		// reset the command size to avoid execution to the serial
		auto cmdsize=0;
	}

	void delete_consumer(const typename consumer_manager_type::key_type & id){
		// Delete a consumer
		// Just the key is needed
		std::cout << progname << " A delete consumer found for key: "<< id << std::endl;
		consumermgr.remove_consumer(id);
		// reset the command size to avoid execution to the serial
		auto cmdsize=0;
	}
	std::string process_msg(const std::string &c){
			int cmdsize=-1;
			std::string msg(c);
			// Create an anonymous sequence to handle sequences not to be recorded
			const std::string ANOMSEQ("noname");
//			seqname=ANOMSEQ;
			// Now
			// Check if is a command for the process a start seqence or wait condition is requested
			// A string transformation is needed to be sure we are using the same case
//			std::transform( msg.begin(), msg.begin()+6, msg.begin(), tolower);
			// Deserialize the process ID and the actions. These are the first two bytes of the message
			id_and_action dec_id_action;
			// Deserialize the process ID and the action, use the automatic deserialization option using a representation deserializer
			auto [lastptr,val]=deserializer_.deserialize(&dec_id_action,msg.c_str());
			// Check if the process_id is the same as the one in the message
			fcutils::print_vals((std::cout<<std::hex<<std::showbase),static_cast<std::uint32_t>(dec_id_action.id),static_cast<std::uint32_t>(dec_id_action.action));
			std::cout<<std::dec<<std::noshowbase;
			if(dec_id_action.id==procid_) {
				std::cout << progname << " right proc_id found !!!" << std::endl;
			}

//			if( ( dec_id_action!=std::string::npos )  ){
//				// a sequence DELETE is found.
//				seqdef_=false;
//				// Extract the sequence name (ID) from the message and remove it
//				// Remember the string first two char are the hex identifier "0x", then there is the address, the command and the ID
//				std::size_t start=2+2+2;
//				seqname=msg.substr(start,proc::cmd::size::SEQIDSIZE*2 );
//				// reset the command size to avoid execution to the serial
//				cmdsize=0;
//				std::cout << progname << " A delete sequence found. Sequence ID: " << seqname << std::endl;
//				// remove the sequence
//				sequencemgr.remove_sequence(seqname);
//			}
//			if( ( msg.find(SEQSTART)!=std::string::npos ) ){
//				  // a sequence is starting just connect to the signal
//				  seqdef=true;
//				  // reset the connection
//				  seqname=ANOMSEQ;
//				  sequencemgr.add_sequence(seqname);
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//				  std::cout << progname << " A start of sequence found" << std::endl;
//			  }
//			  if( ( msg.find(SEQEND)!=std::string::npos )  ){
//				  // a sequence is starting just connect to the signal
//				  seqdef=false;
//				  std::cout << progname << " An end of sequence found. Executing the sequence ..." << std::endl;
//				  // execute the sequence
//				  sequencemgr.execute_sequence(ANOMSEQ);
//				  sequencemgr.remove_sequence(ANOMSEQ);
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if( ( msg.find(SEQSTARTREC)!=std::string::npos ) ){
//				  // Start a sequence recording
//				  seqdef=true;
//				  // Extract the sequence name (ID) from the message and execute it
//				  // Remember the string first two char are the hex identifier "0x", then there is the address, the command and the ID
//				  std::size_t start=2+2+2;
//				  seqname=msg.substr(start,proc::cmd::size::SEQIDSIZE*2 );
//				  // reset the connection
//				  sequencemgr.add_sequence(seqname);
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//				  std::cout << progname << " A start of sequence recording found. Sequence ID: " << seqname << std::endl;
//			  }
//			  if( ( msg.find(SEQENDREC)!=std::string::npos )  ){
//				  // a sequence is starting just connect to the signal
//				  seqdef=false;
//				  std::cout << progname << " An end of sequence recording found. " << std::endl;
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if( ( msg.find(WAIT)!=std::string::npos ) ){
//				  // Check if any duration has been specified
//				  // Create the converter
//				  convtochar< cmdval<proc::cmd::size::SETCMD> > converter;
//				  converter.transform(msg);
//				  cmdval<proc::cmd::size::SETCMD> wc= converter.getcmd();
//				  int d=static_cast<int>(wc[wc.eoc-1]+(wc[wc.eoc-2]<<8));
//				  std::cout << progname << " A wait found for " << d << " milliseconds" << std::endl;
//				  if(seqdef) {
//					  sequencemgr.add_cmd_to_sequence(seqname, boost::bind(&proc_commons::sleep_for,this,d) );
//				  } else {
//					  sleep_for(d);
//				  }
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if( msg.find(SEQEXEC)!=std::string::npos ) {
//				  // execute a sequence
//				  // Extract the sequence name (ID) from the message and execute it
//				  // Remember the string first two char are the hex identifier "0x", then there is the address, the command and the ID
//				  std::size_t start=2+2+2;
//				  sequencemgr.execute_sequence(msg.substr(start,proc::cmd::size::SEQIDSIZE*2 ));
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if(msg.find(ENABLELOG)!=std::string::npos ) {
//				  // Enable a log file of the std output
//				  // Open a log file with the ID and time stamps
//				  logger.enable(std::cout);
//				  std::cout << progname << " An enable log found. Duplicating standard output on log file: "<< logger.get_name() << std::endl;
//				    // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if(msg.find(DISABLELOG)!=std::string::npos ) {
//				  // Disable a log file of the std output
//				  logger.disable(std::cout);
//				  std::cout << progname << " A disable log found. Disable duplication of standard output on log file: "<< logger.get_name() << std::endl;
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if(msg.find(SEQUENCERREC)!=std::string::npos ) {
//				  // Create a new sequencer and store it into the map
//				  // Remember!!! Sizes are in Bytes, we need to multiply by 2 to get the position inside a string.
//				  std::size_t sequenceridpos=1+proc::cmd::size::CMDADDRESSSIZE + proc::cmd::size::CMDACTIONSIZE;
//				  std::size_t seqidpos=sequenceridpos+proc::cmd::size::SEQUENCERIDSIZE;
//				  std::size_t timetype=seqidpos+proc::cmd::size::SEQIDSIZE;
//				  std::size_t timeval=timetype+proc::cmd::size::TIMECNTUNITSIZE;
//				  std::string sequencerid=msg.substr(sequenceridpos*2,proc::cmd::size::SEQUENCERIDSIZE*2);
//				  // Create the converter and transform the string into an array of bytes.
//				  convtochar< cmdval<proc::cmd::size::SEQUENCERCMD> > converter;
//				  converter.transform(msg);
//				  std::string seqid=msg.substr(seqidpos*2,proc::cmd::size::SEQIDSIZE*2);
//				  proc::cmd::byte ttype=converter.getcmd().val()[timetype-1];
//				  // We need to swap the two bytes of the time value, let's do it the hard way
//				  int tval=static_cast<int>(converter.getcmd()[timeval]+(converter.getcmd()[timeval-1]<<8) );
//				  std::cout << progname << " A record a new sequencer found. ID (hex): "
//						  << sequencerid
//						  << ", Sequence to be sequenced ID: " << seqid
//						  << ", Time type: "<< static_cast<unsigned int>(ttype)
//						  << ", Time value: " << tval <<  std::endl;
//				  sequencermgr.add_sequencer(sequencerid,ttype,tval,sequencemgr[seqid]);
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if(msg.find(SEQUENCERSTART)!=std::string::npos ) {
//				  // Start sequencer
//				  // Remember!!! Sizes are in Bytes, we need to multiply by 2 to get the position inside a string.
//				  std::size_t sequenceridpos=1+proc::cmd::size::CMDADDRESSSIZE + proc::cmd::size::CMDACTIONSIZE;
//				  std::string sequencerid=msg.substr(sequenceridpos*2,proc::cmd::size::SEQUENCERIDSIZE*2);
//				  // Create the converter and transform the string into an array of bytes.
//				  std::cout << progname << " A start sequencer found. Sequencer ID (hex): " << sequencerid << std::endl;
//				  if(sequencerid==proc::cmd::ALL){
//					  // Start all
//					  std::cout << progname << " Starting all sequencers" << std::endl;
//					  sequencermgr.start_all_sequencer();
//				  } else {
//					  sequencermgr.start_sequencer(sequencerid);
//				  }
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if(msg.find(SEQUENCERSTOP)!=std::string::npos ) {
//				  // Stop sequencer
//				  // Remember!!! Sizes are in Bytes, we need to multiply by 2 to get the position inside a string.
//				  std::size_t sequenceridpos=1+proc::cmd::size::CMDADDRESSSIZE + proc::cmd::size::CMDACTIONSIZE;
//				  std::string sequencerid=msg.substr(sequenceridpos*2,proc::cmd::size::SEQUENCERIDSIZE*2);
//				  // Create the converter and transform the string into an array of bytes.
//				  std::cout << progname << " A stop sequencer found. Sequencer ID (hex): " << sequencerid << std::endl;
//				  if(sequencerid==proc::cmd::ALL){
//					  // Stop all
//					  std::cout << progname << " Stopping all sequencers" << std::endl;
//					  sequencermgr.stop_all_sequencer();
//				  } else {
//					  sequencermgr.stop_sequencer(sequencerid);
//				  }
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//			  if(msg.find(SEQUENCERDEL)!=std::string::npos ) {
//				  // Delete sequencer
//				  // Remember!!! Sizes are in Bytes, we need to multiply by 2 to get the position inside a string.
//				  std::size_t sequenceridpos=1+proc::cmd::size::CMDADDRESSSIZE + proc::cmd::size::CMDACTIONSIZE;
//				  std::string sequencerid=msg.substr(sequenceridpos*2,proc::cmd::size::SEQUENCERIDSIZE*2);
//				  // Create the converter and transform the string into an array of bytes.
//				  std::cout << progname << " A delete sequencer found. Sequencer ID (hex): " << sequencerid << std::endl;
//				  if(sequencerid==proc::cmd::ALL){
//					  // Stop all
//					  std::cout << progname << " Deleting all sequencers" << std::endl;
//					  // stop all sequences before removing them
//					  sequencermgr.stop_all_sequencer();
//					  sequencermgr.remove_all_sequencer();
//				  } else {
//					  // Stop the sequence before romoving it
//					  sequencermgr.stop_sequencer(sequencerid);
//					  sequencermgr.remove_sequencer(sequencerid);
//				  }
//				  // reset the command size to avoid execution to the serial
//				  cmdsize=0;
//			  }
//
//			  if(msg.find(PROCSPAWNCMD)!=std::string::npos ) {
//				  // Just spawn the command send in the string
//				  // Store the command in a temporary string
//				  std::string s(msg,proc::cmd::position::SPAWNCMDPOS);
////				  //Launch child process
////				  std::cout << progname << " Spawing command: " << s << std::endl;
////				  if(0 != std::system(s.c_str()))
////					  std::cout << progname << " ERROR ->>> Error spawing command: " << s << std::endl;
//				  if(seqdef) {
//					  sequencemgr.add_cmd_to_sequence(seqname, boost::bind(&proc_commons::spawn_cmd,this,s) );
//				  } else {
//					  spawn_cmd(s);
//				  }
//				  // reset the command size to avoid dispatch to the serial port
//				  cmdsize=0;
//			  }
//
//			  if(msg.find(PROCECHOCMD)!=std::string::npos ) {
//				  // Just echo the msg received, as it is, adding the PROC ID
//				  std::cout << progname << " Echo command: " << msg << std::endl;
//				  msg=fcutils::string::add_procid(msg);
//				  // reset the command size to avoid dispatch to the serial port
//				  cmdsize=0;
//			  }
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
//
//			  reply_=msg;
//			  cmdexec=true;
//			  switch(cmdsize){
//				  case 0:
//				  	  {
//				  		  // A sequence start or end, or a wait message was found
//				  		  reply_+=" done \n";
//				  		  // No commands to be executed
//				  		  cmdexec=false;
//				  	  }
//					  break;
//				  case -1:
//				  {
////					  std::cout << progname <<" No PROC commands found" << std::endl;
//					  // The command can be executed
//					  cmdexec=true;
//				  }
//				  break;
//				  default:
//					  std::cout << progname << " Error !!!!!!!!! Wrong command size -->>>> " << cmdsize << std::endl;
//					  cmdexec=false;
//					  break;
//			  }
//			  return reply_;
			return "ok";
		  }


	// Save some  ant strings to be used for checking the message header
//		std::string PROCHEADER;
//		std::string SEQDEL;
//		std::string SEQSTART;
//		std::string SEQEND;
//		std::string SEQSTARTREC;
//		std::string SEQENDREC;
//		std::string SEQEXEC;
//		std::string WAIT;
//		std::string ENABLELOG;
//		std::string DISABLELOG;
//		std::string OPCMDHEADER;
//		std::string SSCMDHEADER;
//		std::string SEQUENCERREC;
//		std::string SEQUENCERSTART;
//		std::string SEQUENCERSTOP;
//		std::string SEQUENCERDEL;
//		std::string PROCSPAWNCMD;
//		std::string PROCECHOCMD;
//		std::string PROCENABLENETLOG;
//		std::string PROCDISABLENETLOG;
//		std::string CONSUMERADD;
//		std::string CONSUMERDEL;
//		std::string reply;

	};

} // namesapce command
} // namespace fcutils


#endif // PROCESS_COMMON_HPP_ 
