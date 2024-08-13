///
///
/// @brief A new file
///
/// @file ostream_logger.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 25 mar 2022 
///

#ifndef OSTREAM_LOGGER_HPP_
#define OSTREAM_LOGGER_HPP_

#include <iosfwd>
#include <chrono>
#include <thread>
#include <string>
#include <ctime>
#include <sstream>
#include <fstream>

#include <boost/iostreams/tee.hpp> // tee_stream & tee_device
#include <boost/iostreams/stream.hpp> // iostream

namespace fcutils{
namespace iostream{

class ostream_logger{
public:
	// Open a log file with the ID and time stamps
	  // Duplicate std::cout using the boost tee stream
	  typedef boost::iostreams::tee_device<std::ostream, std::ofstream> TeeDevice;
	  typedef boost::iostreams::stream<TeeDevice> TeeStream;

	  ostream_logger(std::ostream & o, const char * pre="NoName", const char * post="log", const bool long_name=true ):
		  prefix_(pre), postfix_(post), opencount_(0), long_name_(long_name), out_(NULL) {enable(o);}
	  ostream_logger(const char * pre="NoName", const char * post="log", const bool long_name=true):
		  outputdevice_(NULL), logger_(NULL), prefix_(pre), postfix_(post), opencount_(0), long_name_(long_name), out_(NULL) {};
	  ~ostream_logger(){
//		  reset();
		  if(logger_!=NULL){
			  delete logger_;
		  }
		  if(outputdevice_!=NULL){
			  delete outputdevice_;
		  }
		  if(out_!=NULL){
			  delete out_;
		  }
		  if(fout_.is_open()) {
			  std::cerr << "opps fout_ is still open " << std::endl;
			  fout_.close();
		  }
	  }

	  void enable(std::ostream &o){
		  // Check if is already enabled
		  if( out_==NULL ){
			  std::stringstream ss;
			  ++opencount_;
			  ss << prefix_;
			  // If long name, add information on time and process ID
			  if(long_name_) {
				// Convert the time point into a string YYYYMMDD_HHMM
				std::time_t n = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
				char c[15];
				std::strftime( c, 15, "%Y%m%d", std::localtime(&n) );
				ss<<"_id_"<< std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id())) << '_' <<c;
			  } // if(long_name)
			  ss << "_" << opencount_ << "." << postfix_;
			  fout_name_=ss.str();
			  // Close the opened fstream and delete the tee and the ostream copy
			  reset();
			  fout_.open(fout_name_.c_str());
			  // Temporary copy the std::cout buffer
			  out_=new std::ostream(o.rdbuf());
			  outputdevice_=new TeeDevice(*out_, fout_);
			  logger_=new TeeStream(*outputdevice_);
			  o.rdbuf(logger_->rdbuf());
			  //		  std::cout << PROGNAME << " Duplicating standard output on log file: "<< ss.str() << std::endl;
		  }
	  }
	  void disable(std::ostream &o){
		  if(out_!=NULL) o.rdbuf(out_->rdbuf());
		  reset();
	  }

	  bool is_open() const { return fout_.is_open(); }
	  const char * get_name() const { return (fout_name_.size()>0)?fout_name_.c_str():nullptr;}

private:
	  TeeDevice * outputdevice_;
	  TeeStream * logger_;
	  int opencount_;
	  bool long_name_;
	  std::ofstream fout_;
	  std::ostream * out_;
	  std::string fout_name_;
	  std::string prefix_;
	  std::string postfix_;
	  void reset(){
		  if(logger_!=NULL){
			  logger_->flush();
			  logger_->close();
			  delete logger_;
			  logger_=NULL;
		  }
		  if(outputdevice_!=NULL){
			  delete outputdevice_;
			  outputdevice_=NULL;
		  }
		  if(out_!=NULL){
			  delete out_;
			  out_=NULL;
		  }
		  if(fout_.is_open()) {
			  std::cerr << "opps fout_ is still open " << std::endl;
			  fout_.close();
		  }

	  }
};

} // nsmespace iostream
} // namespace fcutils



#endif // OSTREAM_LOGGER_HPP_ 
