///
///
/// @brief A new file
///
/// @file posix_message_queue.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 14 giu 2021 
///

#ifndef POSIX_MESSAGE_QUEUE_HPP_
#define POSIX_MESSAGE_QUEUE_HPP_

extern "C" {
#include <mqueue.h> // For message queue related definition
#include <signal.h> // for sigevent, sigval_t etc. etc.
};

#include <iostream>
#include <string>
#include <cstring> // for strerror
#include <chrono> // for duration and clock utilities

namespace {
 template < class Int > bool test_bit(const Int f, const Int b){
	 return (f&b)==b;
 }

 //typedef union sigval sigval_t;
 using sigval_t=union sigval;
}

namespace fcutils{
namespace system{

class posix_message_queue{
public:
	enum class ERRORS: int { BAD=-1, GOOD=0 };

	using descriptor_type=mqd_t;
	using attribute_type=mq_attr;
	using permission_type=mode_t;
	using flag_type=int;
	using priority_type=unsigned int;
	using message_type=char;

	posix_message_queue(const char *n="", const flag_type f=flag_type{O_RDWR|O_CREAT|O_NONBLOCK}, const permission_type p=permission_type{S_IRUSR|S_IWUSR}, const attribute_type & a=attribute_type{0L, -1L, -1L, 0L} ):
		name_{n}, flag_{f}, permission_{p}, attribute_{a} {};

	~posix_message_queue() {
		if(opened_) mq_close(descriptor_);
		if(unlink_on_destruction_) mq_unlink(name_.c_str());
	}

	void set_name(const std::string & s){ name_=s;}
	void set_attr(const attribute_type & a) { attribute_=a; }
	void set_perm(const permission_type & p) { permission_=p; }
	void set_unlink(bool b=true){ unlink_on_destruction_=b; }

	const char * get_name() const { return name_.c_str(); }
	const descriptor_type get_descriptor() const { return descriptor_;}
	const attribute_type & get_attribute() const { if(opened_) update_attribute(); return attribute_; }
	const std::size_t get_msgsize() const { return attribute_.mq_msgsize; }
	const flag_type get_flag() const { return flag_;}
	const permission_type get_permission() const { return permission_;}

	bool is_open() const { return opened_; }

	bool close(){
		int i=mq_close(descriptor_);
		bool success=(i==static_cast<int>(ERRORS::GOOD));
		opened_ = (success)? false: true;
		return success;
	}


	bool open(const std::string & name, const flag_type f){
		name_=name;
		flag_=f;
		auto success=open();
		return success;
	};


	bool open(const std::string & name, const flag_type f, const permission_type & p, const attribute_type & a){
		name_=name;
		flag_=f;
		attribute_=a;
		permission_=p;
		auto success=open();
		return success;
	};

	bool open(){
		// Test if queue name has been set, otherwise return false
		if(name_.size()==0) return false;
		// Test if queue creation is requested. If so the permission must be added along well the attributes in the mq_open call
		if( test_bit(flag_,O_CREAT) ) {
			// Test if the attribute message size has been changed otherwise use the system default
			auto ptr=( (attribute_.mq_maxmsg==-1l) && (attribute_.mq_msgsize==-1l) )? nullptr: & attribute_;
			descriptor_ = mq_open(name_.c_str(), flag_,permission_, ptr);
		} else {
			descriptor_=mq_open(name_.c_str(), flag_);
		}
		if (static_cast<int>(descriptor_) == static_cast<int>(ERRORS::BAD)) {
			std::cout << "Error in opening queue: " << name_
						<< strerror(errno) << std::endl;
			return false;
		}
		if(!update_attribute()) return false;
		opened_=true;
		return true;
	}

	bool send(const char * msg, std::size_t length, int priority=0){
		if(!opened_) return false;
		int temp=mq_send(descriptor_, msg, length, priority);
		return (temp==static_cast<int>(ERRORS::GOOD));
	}

	bool send(const std::string & msg, priority_type priority=0){
		return send(msg.c_str(),msg.size(),priority);
	}

	template <class Duration=std::chrono::seconds > bool timedsend(message_type *msg_ptr, std::size_t msg_len,  const Duration & d=Duration{1}, priority_type priority=0){
		if(!opened_) return false;
		if(msg_len<attribute_.mq_msgsize) return false;
		// Calculate the absolute time
		auto now=std::chrono::system_clock::now();
		now=now+d;
		timespec timeout{std::chrono::system_clock::to_time_t(now),0};
		int temp=mq_timedsend(descriptor_, msg_ptr, msg_len, priority, &timeout);
		return (temp==static_cast<int>(ERRORS::GOOD));
	}

	ssize_t receive(message_type *msg_ptr, std::size_t msg_len, priority_type *priority=nullptr){
		if(!opened_) return static_cast<ssize_t>(ERRORS::BAD);
		if(msg_len<attribute_.mq_msgsize) return static_cast<ssize_t>(ERRORS::BAD);
		return mq_receive(descriptor_, msg_ptr, msg_len, priority);
	}

	template <class Duration=std::chrono::seconds > ssize_t timedreceive(message_type *msg_ptr, std::size_t msg_len,  const Duration & d=Duration{1}, priority_type *priority=nullptr){
		if(!opened_) return static_cast<ssize_t>(ERRORS::BAD);
		if(msg_len<attribute_.mq_msgsize) return static_cast<ssize_t>(ERRORS::BAD);
		// Calculate the absolute time
		auto now=std::chrono::system_clock::now();
		now=now+d;
		timespec timeout{std::chrono::system_clock::to_time_t(now),0};
		return mq_timedreceive(descriptor_, msg_ptr, msg_len, priority, &timeout);
	}

	/// This notify can be used to register the posix_message_queue object in the callback data (in sigval_t)
	/// In this way the function can re-notify, provided that the object is still alive ...
	/// The queue is set in THREAD notification mode.
	bool notify(void (*f)(sigval_t))
	{
		if(!opened_) return false;
		sigevent sev;
		sev.sigev_notify = SIGEV_THREAD; // Notify via thread
		sev.sigev_notify_function = f;
		sev.sigev_notify_attributes = NULL;
	// Could be pointer to pthread_attr_t structure
		sev.sigev_value.sival_ptr = this; // Argument to threadFunc()
		auto temp=mq_notify(descriptor_, &sev);
		if (temp==static_cast<int>(ERRORS::BAD)){
			std::cout << "Error in mq_notify: "
						<< strerror(errno) << std::endl;
			return false;
		}
		return true;
	}

	/// This notify gives you masimum flexibility in specification of the sigevent object.
	bool notify(const sigevent & sev)
	{
		if(!opened_) return false;
		auto temp=mq_notify(descriptor_, &sev);
		if (temp==static_cast<int>(ERRORS::BAD)){
			std::cout << "Error in mq_notify: "
						<< strerror(errno) << std::endl;
			return false;
		}
		return true;
	}

private:
	 std::string name_;
	 descriptor_type descriptor_;
	 mutable attribute_type attribute_;
	 permission_type permission_;
	 flag_type flag_;
	 bool unlink_on_destruction_{false};
	 bool opened_{false};
	 bool update_attribute() const {
		auto temp= mq_getattr(descriptor_, &attribute_);
		if ( temp == static_cast<int>(ERRORS::BAD)) {
			std::cout << "Error in getting queue: " << name_<< " attributes. Error: "
					<< strerror(errno) << std::endl;
			return false;
		}
		return true;
	 }
};

} // namespace system
} // namespace fcutils



#endif // POSIX_MESSAGE_QUEUE_HPP_ 
