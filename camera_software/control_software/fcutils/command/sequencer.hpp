///
///
/// @brief A sequencer class. It execute a given action according to a given frequency or a given number of time.
///
/// @file sequencer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 17 mar 2022 
///

#ifndef SEQUENCER_HPP_
#define SEQUENCER_HPP_

#include <iostream>
#include <string>
#include <chrono>
#include <atomic>
#include <memory>
#include <thread>
#include <functional>
#include <variant>

#include "fcutils/utility/overload.hpp"
#include "fcutils/concurrency/guarded_thread.hpp"
#include "fcutils/chrono/chrono.hpp"

namespace {
	using namespace std::chrono_literals;
}

namespace fcutils{
namespace command{

class base_sequencer{
public:
	using count_type		= std::size_t;
	base_sequencer(){};
	base_sequencer( const base_sequencer &) = delete;
	base_sequencer & operator=(const base_sequencer &) = delete;
	base_sequencer( base_sequencer &&) noexcept = default;
	base_sequencer & operator=( base_sequencer &&) noexcept = default;
	virtual ~base_sequencer(){};

	virtual bool is_started() const = 0;
	virtual bool is_time() const =0;
	virtual bool is_count() const = 0;
	virtual void set_interval( const int & ) = 0;
	virtual void set_repetition( const count_type &  ) = 0;
	virtual void start()= 0;
	virtual void stop()=0;
	virtual void print( std::ostream & o) const = 0;
};

// To execute a given list of sequences with a given time interval or a given number of repetition
// Is user responsibility to handle the action lifetime
// the default behavior is a time sequencer, i.e. it repeats the action and waits until is stopped.
// If the IsTime is set to false, it consider the Duration interval as a counter and repeat the action for the counter number of time or until is stopped

template<typename Action = std::function<void () >, typename Duration =std::chrono::milliseconds>
class sequencer: public base_sequencer  {
public:
	using base_type			= base_sequencer;
	using self_type			= sequencer <Action, Duration>;
	using action_type 		= Action;
	using action_ptr_type	= std::shared_ptr<action_type>;
	using duration_type 	= Duration;
	using count_type		= base_sequencer::count_type;
	using data_type			= std::variant<duration_type, count_type>;
	using thread_type		= fcutils::concurrency::joinable_guarded_thread;

	// A creator with a duration. This is the default creator
	template <typename D> sequencer( const D & d=30s, action_ptr_type action=NULL ): started_{false}, action_{action} {
		// Check if D is an integer type, if not, assume is a chrono::duration type
		if constexpr (std::is_integral<D>::value){
			data_=static_cast<count_type>(d);
		} else {
			data_=std::chrono::duration_cast<duration_type>(d);
		}
	};

	sequencer(const sequencer & s) = default;
	sequencer & operator=(const sequencer & s) = default;
	sequencer(sequencer && s) noexcept = default;
	sequencer & operator=(sequencer && s) noexcept = default;

	~sequencer() {
		std::cout << "Oops destroying sequencer of type: " << std::string{(is_time())?"":"duration"} << std::string{(is_count())?"":"count"} << std::endl;
		stop();
	}

	bool is_time() const override {
		return std::holds_alternative<duration_type>(data_);
	}
	bool is_count() const override {
		return std::holds_alternative<count_type>(data_);
	}
	bool is_started() const override {
		return started_.load();
	}

	void start() override {
		if (!action_)
			return;
		if (!started_.load()) {
			started_ = true;
			// It should never be joinable at this point !!!
			if (sequencer_thread_.joinable())
				sequencer_thread_.join();
//			auto ovl= [this]( auto & v ) mutable{ execute(v);};
			// Launch the thread ...
//			std::visit(ovl,data_);
			sequencer_thread_ = std::thread{
					[this](){
						std::visit(
								[this](auto & v) mutable { execute(v);}
								,data_); // std::visit
						} // Lambda for thread
			}; // std::thread

//			sequencer_thread_ = std::thread([this](){std::visit(ovl_,data_);});
//			sequencer_thread_ = boost::thread(
//					boost::bind(
//							static_cast<void (self_type::*)(Int2Type<IsTime>)> (&self_type::execute),
//							this,
//							Int2Type<IsTime>() )
//							);
		}
	}
	void stop() override {
		if(!action_) return;
		if(started_.load()) {
			started_= false;
			// interrupt the thread and wait for it
//			sequencer_thread_.interrupt();
			sequencer_thread_.join();
		}
	}


	// Set an interval
	template <typename D = duration_type>
	void set_interval( const D & rel_time ){
		// restart the sequencer after changing the interval, if the sequencer was running
		bool temp=started_.load();
		if(started_.load()) stop();
		data_= std::chrono::duration_cast<duration_type>(rel_time);
		if(temp) start();
	}

	void set_interval( const int & rel_time ) override {
		set_interval(duration_type{rel_time});
	}

	void set_repetition( const count_type & c ) override {
		// restart the sequencer after changing the interval, if the sequencer was running
		bool temp=started_.load();
		if(started_.load()) stop();
		data_= c;
		if(temp) start();
	}

	void set_action(action_ptr_type action) {
		action_=action;
		if(started_) {
			stop();
			start();
		}
	}

	void print( std::ostream & o =std::cout ) const override {
		o << " Sequencer of type " << (is_time()?"time: ":"") << (is_count()?"count: ":"");
		auto ovl=fcutils::overload{
			[&o]( const count_type & c ) { o << c << '\n';},
			[&o]( const duration_type & d){ o << d.count() << fcutils::chrono::suffix<duration_type>() << '\n'; }
		};
		std::visit( ovl, data_);
	}

private:
	data_type data_;
	std::atomic<bool> started_{false};
	action_ptr_type action_;
	thread_type sequencer_thread_;

	void execute(duration_type & d) {
		if(d>std::chrono::duration_cast<duration_type>(0ms)){
			while(started_.load()) {
				action_->operator()();
				std::this_thread::sleep_for(d);
			}
		}
	}

	void execute(count_type & d) {
		if(!started_.load()) return;
		auto temp=d;
		while(started_.load() && (temp>0) ) {
			action_->operator()();
			--temp;
		}
		started_=false;
	}
};

} // namespace command
} // namespace fcutils


#endif // SEQUENCER_HPP_ 
