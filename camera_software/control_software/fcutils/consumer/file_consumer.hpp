///
///
/// @brief A new file
///
/// @file file_consumer.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 14 mar 2022 
///

#ifndef FILE_CONSUMER_HPP_
#define FILE_CONSUMER_HPP_

#include <fstream>
#include <iterator>
#include <chrono>
#include <ctime>
#include <thread>
#include <string>

#include "fcutils/checker/true_checker.hpp"
#include "fcutils/consumer/consumer_base.hpp"

namespace{
	using namespace std::chrono_literals;
}


namespace fcutils{
namespace consumer{

// The construction string contains the file type and the option to logrotate or not.
// By default it rotates the log every hour.
// if % is present after the filename, it will log rotate for the amount, in mininutes, after the % char.
// Examples:
// filename%100, it will create the filename, in binary mode and log rotate it every 100minutes.

enum class FILE_FORMAT: bool { BINARY=true, TEXT=false };
template <typename Data=unsigned char, template <typename D> class Checker = fcutils::checker::true_checker, FILE_FORMAT Format=FILE_FORMAT::TEXT,
		char EOL='\0', typename Duration=std::chrono::minutes>
class file_consumer: public consumer_base<Data,Checker> {

//class binaryfileconsumer: public consumer_base{
public:
	using self_type			=file_consumer<Data,Checker,Format,EOL,Duration>;
	using consumer_base_type=consumer_base<Data,Checker>;
	using checker_type 		=typename consumer_base_type::checker_type;
	using data_type			=typename consumer_base_type::data_type;
	using file_type			=std::ofstream;
	using duration_type		=Duration;
	using time_point_type	=std::chrono::time_point<std::chrono::system_clock>;
	using clock_type 		=time_point_type::clock;

	template<typename C=checker_type>
	file_consumer(const char * name, std::size_t timeout=15, bool longname=true, const char * ext="dat", C && checker=checker_type{}):
	file_consumer(name,duration_type{timeout},longname,ext,checker){};
//	consumer_base_type{std::forward<checker_type>(checker)}, logrotate_{(duration_type{timeout}>0min)?true:false}, timeout_{duration_type{timeout}},opencount_{0},
//	longname_{longname}, prefix_{name}, postfix_{ext}, start_{}
//	{
//		// Just calculate the hash of the thread id, to be used in case of long name. Please note that is going to be the thread of the creation block. It could be different from the execution thread.
//		thread_id_=std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id()));
//	};

	template<typename C=checker_type>
	file_consumer(const char * name, duration_type timeout=15min, bool longname=true, const char * ext="dat", C && checker=checker_type{}):
		consumer_base_type{std::forward<checker_type>(checker)}, logrotate_{(timeout>0min)?true:false}, timeout_{timeout},opencount_{0},
		longname_{longname}, prefix_{name}, postfix_{ext}, start_{}
		{
			// Just calculate the hash of the thread id, to be used in case of long name. Please note that is going to be the thread of the creation block. It could be different from the execution thread.
			thread_id_=std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id()));
		};
	//	Move constructors
	file_consumer( file_consumer && ) noexcept = default;
	self_type & operator=( file_consumer && ) noexcept = default;

	~file_consumer(){
		if(fout_.is_open()) {
			close();
		}
	}


	self_type & operator()(const data_type * c, std::size_t length) override {
		// Create the socket, open it and bind it to the local endpoint
		if(is_expired()) open();
		if constexpr (Format == FILE_FORMAT::TEXT){
			std::ostreambuf_iterator<data_type> oo (fout_.rdbuf()); // output iterator for cout
			std::copy(c,c+length,oo);
		} else {
			fout_.write(c,length);
		}
		if constexpr(EOL!='\0') {
			auto temp=EOL;
			fout_.write(&temp,sizeof(decltype(EOL)));
		}
		fout_.flush();
		consumer_base_type::update_stats(length);
		return *this;
	}

	void reset() override {
		opencount_=0;
		consumer_base_type::reset();
	}
	bool is_logrotate() const { return logrotate_;}
	bool is_longname() const { return longname_;}
	bool is_binary() const { return static_cast<bool>(Format); }

	// Set a default log rotation of 15 minutes
	template <typename LogDuration = duration_type>
	void set_logrotation( const LogDuration & d=std::chrono::duration_cast<duration_type>(15min)){
		logrotate_=true;
		timeout_=std::chrono::duration_cast<duration_type>(d);
		set_start();
	}
	// Set a log duration in minutes not using the chrono duration. Relies on implicit casting
	void set_logrotation( int d=15 ){
		logrotate_=true;
		timeout_=duration_type{d};
		set_start();
	}

	void stop_logrotation() { logrotate_=false; }

	void set_longname() { longname_=true; }
	void stop_longname() { longname_=false; }

	const char * get_name() const { return fout_name_.c_str();}

private:
	std::string fout_name_;
	std::ofstream fout_;
	bool logrotate_{true};
	duration_type timeout_;
	time_point_type start_;
	int opencount_;
	bool longname_;
	std::string prefix_;
	std::string postfix_;
	std::string thread_id_;


	bool is_expired ()
	{
		// if logrotate is not enabled return the file status
		if ( ! logrotate_) return true;
		// If logrotate is enable, then check for elapsed time and let the open method do the job of checking the status
		duration_type elapsed=std::chrono::duration_cast<duration_type>(clock_type::now()-start_);
		return (elapsed> timeout_);
	}

	void set_start(){
		start_= clock_type::now();
	}

	void open(){
		// Close the opened fstream if any
		close();
		// Create the string with the new name initializing it with the prefix
		std::string ss{prefix_};
		// Increment the counter of open operations
		++opencount_;
		  // Check if logrotate is enabled
		if(logrotate_) {
			// If long name, add information on time and process ID
			if(longname_) {
				// Convert the time point into a string YYYYMMDD_HHMM
				std::time_t n = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
				char c[15];
//				std::strftime( c, 15, "%Y%m%d_%H%M", std::localtime(&n) );
				std::strftime( c, 15, "%Y%m%d", std::localtime(&n) );
				ss+="_id_";
				// Get the id, has to be "hashed" to an integer
				ss+=thread_id_;
				ss+='_';
				// Add the time in minutes
				ss+=c;
			}
			ss+='_';
			ss+=std::to_string(opencount_);
			set_start();
		}
		ss+='.';
		ss+=postfix_;
		fout_name_=ss;
		if constexpr ( Format == FILE_FORMAT::BINARY ) {
//			std::cout << " Opening binay file: " << fout_name_ << std::endl;
			fout_.open(fout_name_.c_str(),std::fstream::binary);
		} else {
//			std::cout << " Opening txt file: " << fout_name_ << std::endl;
			fout_.open(fout_name_.c_str());
		}
	}

	void close(){
		 if(fout_.is_open()) {
			 fout_.close();
		 }
	 }
};

template<typename Data=unsigned char, template <typename D> class Checker=fcutils::checker::true_checker>
using binaryfile_consumer = file_consumer<Data,Checker,FILE_FORMAT::BINARY>;
template<typename Data=unsigned char, template <typename D> class Checker=fcutils::checker::true_checker, char EOL='\n'>
using textfile_consumer = file_consumer<Data,Checker,FILE_FORMAT::TEXT,EOL>;

} // namesapce consumer
} // namesapce fcutils



#endif // FILE_CONSUMER_HPP_ 
