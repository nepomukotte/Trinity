///
///
/// @brief A new file
///
/// @file sequencermanager.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 18 mar 2022 
///

#ifndef SEQUENCERMANAGER_HPP_
#define SEQUENCERMANAGER_HPP_
/*
 * sequencermanager.hpp
 *
 *  Created on: 09/apr/2014
 *      Author: cafagna
 */

#include <map>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <chrono>
#include <type_traits>

#include "fcutils/concurrency/lockable_container.hpp"

//#include <boost/thread/mutex.hpp>
//#include <boost/shared_ptr.hpp>


namespace {
  using namespace std::chrono_literals;
}

namespace fcutils{
namespace command{

template <class T> void deleter(T mp) {
	delete mp.second;
}

const std::string SMPROGNAME(" SEQENCERMANAGER -->> ");


// This is the manager of the base_sequencer classes.
template <template <typename> class Sequencer, typename BaseSequencer, typename ID=std::string>
class sequencermanager{
public:
	using id_type 				= ID;
	template <typename Action> using sequencer_type = Sequencer<Action>;
	using sequencer_base_type 	= BaseSequencer;
//	using sequencer_base_type	= typename
//			std::conditional<
//				std::is_polymorphic<sequencer_type>::value, // Check if is a derived class
//				decltype(typename sequencer_type::base_type, sequencer_type), // if is a derived class, check if a base_type is defined, otherwise use the sequencer_type
//				sequencer_type>::type; // Use the sequencer_type, if is not a derived class
	using sequencer_ptr_type 	= std::shared_ptr<sequencer_base_type>;
	using sequencer_map_type 	= fcutils::concurrency::lockable_container<std::map,id_type,sequencer_ptr_type>;

	// Add a sequencer with identifier id, duration val and action
	template<class Duration, class Action > bool add_sequencer(const id_type & id, const Duration &val, std::shared_ptr<Action> a  ){
		seqmap_[id]=sequencer_ptr_type(new sequencer_type<Action>(val,a));
		return true;
	}

	// Add a sequencer with identifier id
//	template< class Action > bool add_sequencer(const T & id, const proc::cmd::byte duration_type, const duration_val_type &val, Action *a  ){
//	template< class Action > bool add_sequencer(const T & id, const proc::cmd::byte duration_type, const duration_val_type &val, boost::shared_ptr<Action> a  ){
//		switch(duration_type){
//		case proc::cmd::units::COUNTS:
//			seqmap_[id]=create_sequencer<Action>(procu::counts(val),a);
//			break;
//		case proc::cmd::units::MSEC:
//			seqmap_[id]=create_sequencer<Action>(procu::milliseconds(val),a);
//			break;
//		case proc::cmd::units::SEC:
//			seqmap_[id]=create_sequencer<Action>(procu::seconds(val),a);
//			break;
//		case proc::cmd::units::MIN:
//			seqmap_[id]=create_sequencer<Action>(procu::minutes(val),a);
//			break;
//		case proc::cmd::units::HOUR:
//			seqmap_[id]=create_sequencer<Action>(procu::hours(val),a);
//			break;
//		default:
//			std::cout << PROGNAME << " Error !!!! Units is not defined for sequencer :" << id << ", Unit: " << duration_type << std::endl;
//			break;
//		}
//	}

	bool add_sequencer(const id_type &id, sequencer_ptr_type s) {
		bool temp=false;
		seqmap_[id]=s;
		temp=true;
		return temp;
	}

// remove all sequncer
	bool remove_all_sequencer(){
		if( !seqmap_.empty() ) {
			// Found it delete the pointer and erase it
//			deleter(*i);
			seqmap_.clear();
			return true;
		}
		return false;
	}

	// Remove a sequencer with identifier id
	bool remove_sequencer(const id_type & id){
		// Search for it and then delete
		auto i=seqmap_.find(id);
		if( i!=seqmap_.end() ) {
			// Found it delete the pointer and erase it
//			deleter(*i);
			seqmap_.erase(i);
			return true;
		}
		return false;
	}
	// Execute a given sequence
	bool start_sequencer(const id_type & id) {
		// Search for it and start
		auto i=seqmap_.find(id);
		if(i!=seqmap_.end()) {
			// Found it!!! Start the sequencer
			std::cout << SMPROGNAME << " Starting the sequencer :" << id << std::endl;
			(i->second)->start();
			return true;
		}
		return false;
	}

	bool stop_sequencer(const id_type & id) {
			// Search for it and stop
			auto i=seqmap_.find(id);
			if(i!=seqmap_.end()) {
				// Found it!!! Stop the sequencer
				std::cout << SMPROGNAME << " Stopping the sequencer :" << id << std::endl;
				(i->second)->stop();
				return true;
			}
			return false;
	}

	bool start_all_sequencer(){
		bool executed = false;
		for(auto & [key,ptr]: seqmap_){
			ptr->start();
		}
		executed=true;
		return executed;
	}

	bool stop_all_sequencer(){
		bool executed=false;
		for(auto & [key,ptr]: seqmap_){
			ptr->stop();
		}
		executed=true;
		return executed;
	}

	std::size_t size() const { return seqmap_.size(); }

	typename sequencer_map_type::mapped_type & operator[]( const id_type & id){ return seqmap_[id]; }
	const typename sequencer_map_type::mapped_type & operator[]( const id_type & id) const { return seqmap_[id]; }

	void print(std::ostream & o=std::cout ) const {
		o << SMPROGNAME << " Stored sequencers: " << size() << '\n';
		for( const auto & [key,ptr]: seqmap_){
			ptr->print(o);
		}
		o << std::endl;
	}

private:
	sequencer_map_type seqmap_;

};

} // namespace sequencer
} // namespace proc





#endif // SEQUENCERMANAGER_HPP_ 
