///
///
/// @brief A sequence class
///
/// @file sequence.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 16 mar 2022 
///

#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include <iostream>
#include <boost/signals2.hpp>

namespace fcutils{
namespace command{

template <typename IdType=std::size_t, typename Signal = boost::signals2::signal<void ()>, bool Debug=false >
class sequence {
public:
	using self_type 	= sequence<IdType,Signal,Debug>;
	using signal_type 	= Signal;
	using id_type		= IdType;
	using count_type	= std::size_t;
	using command_type	= typename signal_type::slot_type;

	sequence(const id_type & seqid): seqid_{seqid} {};

	sequence(const self_type & ) = default;
	self_type & operator=( const self_type & ) = default;

	sequence( self_type && ) noexcept = default;
	self_type & operator=( self_type &&) noexcept = default;

	~sequence(){
		if constexpr (Debug){
			std::cout << seqid_ << " Killed !!! " << std::endl;
		}
	}

	void addcmd(const command_type & cmd){ seq_.connect(cmd);}
	// execute the sequence incrementing the sequence count
	void execute(){ ++seqcount_; seq_(); }
	// execute a sequence if has been called less than test times
	bool execute_until(count_type test){
		bool done=(seqcount_>=test);
		if(!done) execute();
		return done;
	};

	// Execute if a test is true
	bool execute_it(bool test){
		if(test) execute();
		return test;
	}

	void print(std::ostream & o=std::cout) {
		o << " Sequence ID: " << seqid_ << ", executed: " << seqcount_ << " times, contains: " << seq_.num_slots() << " commands." << std::endl;
	}

private:
	id_type seqid_;
	count_type seqcount_{0};
	signal_type seq_;
};

} // namespace command
} // namespace fcutils



#endif // SEQUENCE_HPP_ 
