///
///
/// @brief A new file
///
/// @file net_duplicator.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 24 mar 2022 
///

#ifndef NET_DUPLICATOR_HPP_
#define NET_DUPLICATOR_HPP_
#include <iosfwd>

#include <boost/iostreams/tee.hpp> // tee_stream & tee_device
#include <boost/iostreams/stream.hpp> // iostream

#include "fcutils/iostream/net_sink.hpp"

namespace fcutils{
namespace iostream{

class base_net_duplicator{
public:
	base_net_duplicator(){};
	virtual ~base_net_duplicator(){};
	base_net_duplicator( const base_net_duplicator &) = delete;
	base_net_duplicator & operator=(const base_net_duplicator &) = delete;
	virtual void enable(std::ostream &o) =0;
	virtual void disable(std::ostream &o) =0;
	virtual const char * get_address() const =0;
	virtual int get_port() const =0;
	virtual bool is_udp() const =0;
	virtual bool is_tcp() const =0;
};

template<class NetSink> class ostream_net_duplicator: public base_net_duplicator{
public:
	  // Duplicate std::ostream on an UDP port using the boost tee stream and ASIO
	using sink_type	= NetSink;
	using TeeDevice	= boost::iostreams::tee_device<std::ostream, sink_type>;
	using TeeStream	= boost::iostreams::stream<TeeDevice>;

	ostream_net_duplicator(std::ostream & o, const char * address="127.0.0.1", const int port=514 ):
		  sink_(address, port), opencount_(0) {enable(o);}

	ostream_net_duplicator(const char * address="127.0.0.1", const int port=514 ):
		sink_(address, port), opencount_(0), out_(NULL), outputdevice_(NULL),duplicator_(NULL) {};

	~ostream_net_duplicator(){
//		reset();
		if(duplicator_!=nullptr){
			delete duplicator_;
		}
		if(outputdevice_!=nullptr){
			delete outputdevice_;
		}
		if(out_!=nullptr){
			delete out_;
		}
	 }
	void enable(std::ostream &o) override {
		++opencount_;
		// Close the opened fstream and delete the tee and the ostream copy disabling it first
//		reset();
		// delete() takes care of resurce cleanup.
		disable(o);
		// Temporary copy the std::cout buffer
		out_=new std::ostream(o.rdbuf());
		outputdevice_=new TeeDevice(*out_, sink_);
		duplicator_=new TeeStream(*outputdevice_);

		o.rdbuf(duplicator_->rdbuf());
	}

	void disable(std::ostream &o) override {
		if(out_!=nullptr) o.rdbuf(out_->rdbuf());
		reset();
	 }
	const char * get_address() const override { return sink_.get_address();}
	int get_port() const override { return sink_.get_port();}
	bool is_udp() const override { return sink_.is_udp; }
	bool is_tcp() const override { return sink_.is_tcp; }
private:
	TeeDevice * outputdevice_;
	TeeStream * duplicator_;
	sink_type sink_;
	std::ostream * out_;
	int opencount_;
	void reset(){
		if(duplicator_!=nullptr){
			if(duplicator_->operator->()){
				duplicator_->flush();
				duplicator_->close();
			}
			delete duplicator_;
			duplicator_=nullptr;
		}
		if(outputdevice_!=nullptr){
			delete outputdevice_;
			outputdevice_=nullptr;
		}
		if(out_!=nullptr){
			delete out_;
			out_=nullptr;
		}
	 }

};

using ostream_udp_duplicator=ostream_net_duplicator<fcutils::iostream::udp_sink> ;
using ostream_tcp_duplicator=ostream_net_duplicator<fcutils::iostream::tcp_sink> ;

} //namespace iostream
} // namespace fcutils



#endif // NET_DUPLICATOR_HPP_ 
