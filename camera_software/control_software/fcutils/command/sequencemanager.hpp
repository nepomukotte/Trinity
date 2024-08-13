///
///
/// @brief A sequence manager
///
/// @file sequencermanager.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 16 mar 2022 
///

#ifndef SEQUENCEMANAGER_HPP_
#define SEQUENCEMANAGER_HPP_

#include <map>
#include <algorithm>
#include <iostream>
#include <string>
#include <memory>
#include <type_traits>

//#include "fcutils/command/sequence.hpp"
#include "fcutils/concurrency/lockable_container.hpp"

//namespace {
//// This function just deletes the second object of a pair, supposing is a pointer type. No typechek is done.
//template <class Pair> void deleter(Pair *p) {
//	if constexpr ( std::is_pointer< typename Pair::second_type >::value) {
//		delete p->second;
//	} else {
//		std::cout << " Ooopsss, is not a pointer!" << std::endl;
//	}
//}
//
//} // namespace

namespace fcutils{
namespace command{


const std::string SPROGNAME(" SEQMANAGER -->> ");

// #define SPROGNAME " SEQMANAGER -->> "

template <typename Sequence, bool Debug = false > class seqmanager{
public:
	using self_type			= seqmanager<Sequence,Debug>;
	using sequence_type		= Sequence;
//	using sequence_type		= fcutils::command::sequence<id_type,Debug>;
	using id_type	 		= typename sequence_type::id_type;
	using command_type		= typename sequence_type::command_type;
	using sequence_ptr_type	= std::shared_ptr<sequence_type>;
	using seqmap_type		= fcutils::concurrency::lockable_container<std::map,id_type, sequence_ptr_type>;

	seqmanager() {};
	seqmanager & operator=(const seqmanager &) = default;

	seqmanager( seqmanager && ) noexcept = default;
	seqmanager & operator= (seqmanager && ) noexcept = default;

	// Remember that the manager owns the sequences. It must delete them. Let's relay on the sharable_pointer for that.
	~seqmanager(){};
	// Add a new sequence with identifier id
	void add_sequence(const id_type & id){
		// remember to add a check if an id already exists
		// Please note that sharable_pointer will take care of deleting the pointer to the pair as soon all the references has been released.
		seqmap_[id]=sequence_ptr_type(new sequence_type(id));
	}
	// Remove a sequence with identifier id
	bool remove_sequence(const id_type & id){
		// Search for id and then delete
		auto i=seqmap_.find(id);
		if(i!=seqmap_.end()){
			// Found it delete the pointer and erase it
			seqmap_.erase(i);
			return true;
		}
		return false;
	}
	// Check if a sequence with a given ID is present
	bool is_stored(const id_type & id) const {
		auto i=seqmap_.find(id);
		return  i!=seqmap_.end();
	}

	// Add a new command to the sequence id
	bool add_cmd_to_sequence(const id_type & id, const command_type & cmd) {
		// remember to add a check on the id existence.
		auto i=seqmap_.find(id);
		if(i!=seqmap_.end()){
			(i->second)->addcmd(cmd);
			return true;
		}
		return false;
	}

	// Execute a given sequence
	bool execute_sequence(const id_type & id) {
		// Search for it and execute if present
		auto i=seqmap_.find(id);
		if(i!=seqmap_.end()) {
			// Found it execute
			if constexpr (Debug){
				std::cout << SPROGNAME << " Executing the sequence :" << id << std::endl;
			}
			(i->second)->execute();
			return true;
		}
		return false;
	}
	// Execute a given sequence without checking
	self_type & operator()( const id_type & id) {
		bool executed= execute_sequence(id);
		if(! executed){
			// Found it reset the signal
			std::cout << SPROGNAME << " ERROR!!! Executing the sequence :" << id << std::endl;
		}
		return *this;
	}

//	void execute_sequence_nocheck(const id_type & id) {
//		bool executed=execute_sequence(id);
//		if(!executed) {
//			// Found it reset the signal
//			std::cout << SPROGNAME << " ERROR!!! Executing the sequence :" << id << std::endl;
//		}
//	}

	std::size_t size(){ return seqmap_.size(); }
	auto find(const id_type & id) {
		return seqmap_.find(id);
	}

	const auto find(const id_type & id) const {
		return seqmap_.find(id);
	}

	auto & operator[]( const id_type & id){ return seqmap_[id]; }
	const auto & operator[]( const id_type & id) const { return seqmap_[id]; }
	auto & at(const id_type & id) { return seqmap_.at(id); }
	const auto & at(const id_type & id) const { return seqmap_.at(id); }

	void print(std::ostream & o=std::cout ) const {
		for( const auto & [id,s]: seqmap_ ){
			s->print();
		}
	}

private:
	seqmap_type seqmap_;

};

} // namespace command
} // namespace fcutils




#endif // SEQUENCEMANAGER_HPP_ 
