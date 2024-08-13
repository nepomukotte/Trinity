///
///
/// @brief A test for the sequencer and sequencermanager classes.
///
/// @file test_sequencermanager.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 18 mar 2022 
///


#include <iostream>
#include <memory>

#include "fcutils/command/sequencer.hpp"
#include "fcutils/command/sequencermanager.hpp"

struct executor{
	explicit executor( const std::string & s): cmd{s}{};
	executor & operator=(const std::string & s) { cmd=s; return *this;}
	void operator()(){
		std::cout <<" THis is the command: "<< cmd << std::endl;
	}
	std::string cmd;
};

template<typename Action> using sequencer_type=fcutils::command::sequencer<Action>;
using sequencer_base_type = fcutils::command::base_sequencer;
using sequencermanager_type = fcutils::command::sequencermanager<sequencer_type, sequencer_base_type>;

int main(){
	std::string temp{"nothing yet"};
	using action_ptr_type = std::shared_ptr<executor>;
	auto actiond= std::make_shared<executor>(temp);
	auto actionc= std::make_shared<executor>(temp);
	*actionc="count: " + temp;
	sequencermanager_type sm;
	sm.add_sequencer("d10s",10s, actiond);
	sm.add_sequencer("c10",10, actionc);
	sm.start_all_sequencer();

//	fcutils::command::sequencer<executor> seqd{10s,actiond}, seqc{10,actionc};
//	seqd.start();
//	seqc.start();

	std::cout << " Enter anything: \n";
	while( std::getline(std::cin,temp) ){
		*actiond=temp;
		*actionc="count: "+temp;
		sm.start_all_sequencer();
	}
	sm.print();
	std::cout << " Is over, I am going to stop them ......" << std::endl;
	sm.stop_all_sequencer();
//	seqd.stop();
//	seqc.stop();

	return 0;
}


