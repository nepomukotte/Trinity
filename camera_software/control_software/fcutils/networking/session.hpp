///
///
/// @brief A session example in BOOST ASIO. Adapted from the ASIO examples.
///
/// @file session.h
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 13 ott 2021 
///

#ifndef SESSION_HPP_
#define SESSION_HPP_
#include <iostream>
#include <memory> // for std::enable_shared_from_this
#include <utility> // for std::move
#include <algorithm> // for std::copy
#include <type_traits> // for std::void_t
#include <atomic> // for std::atomic_bool
#include <cstddef> // for std::byte

#include <boost/asio.hpp>

namespace fcutils{
namespace networking{

namespace ba=boost::asio;
namespace baip=ba::ip;

///
/// A session is needed by a server or client, to implement the asynchronous reading or writing from or to a socket.
/// A session embed the actual data ... at the moment
///
// Few consumer examples. See fcutils::consumer for more examples
template<class Data> struct do_nothing{
	void operator()(const Data *, std::size_t){};
};

template<class Data> struct print{
	void consume(const Data* d, std::size_t l){
		operator()(d,l);
	}
	void operator()(const Data* d, std::size_t l){
		std::cout << " Bytes received: " << l <<", values:\n";
		std::ostreambuf_iterator<Data> oo (std::cout.rdbuf()); // output iterator for cout
		std::copy(d,d+l,oo);
		std::cout << std::endl;
	}
};

template<class Data> struct print_byte{
	void consume(const Data* d, std::size_t l){
		operator()(d,l);
	}
	void operator()(const Data* d, std::size_t l){
		std::cout << " Bytes received: " << l <<", values:\n";
		for(int i=0; i<l; ++i) std::cout << std::hex<< std::to_integer<short>(static_cast<std::byte>(static_cast<unsigned char>(*(d+i)))) <<' ';
		std::cout << std::dec << std::endl;
	}
};

// Basic checker
// Here we suppose Data is char ....
template<class Data> struct check_eol{
	void consume(const Data* d, std::size_t l){
		operator()(d,l);
	}
	size_t operator()(const Data* d, size_t bytes) {
		// We read one byte at the time until we get '\n'
		// an alternative way checking the full received array
		//		bool found = ( std::find(d, d + bytes, '\n') < (d + bytes) );
		bool found=(*(d+bytes-1)=='\n');
		return found ? 0 : 1;
	}
};

using tcp=baip::tcp;
using udp=baip::udp;

template <class Protocol=tcp, bool Echo=true, template<typename> class Consumer=print, template<typename> class Checker=std::void_t > class session
  : public std::enable_shared_from_this< session<Protocol,Echo,Consumer,Checker> >
{
public:
	using protocol_type=Protocol;
	using socket_type=typename Protocol::socket;
	using buffer_type=char;
	using consumer_type=Consumer<buffer_type>;
	// Use a checker if not void. If a checker is request, then checker_type is used, otherwise a do_nothing class.
	using checker_template_type=Checker<buffer_type>;
	using checker_type= typename std::conditional<std::is_void< checker_template_type >::value,do_nothing<buffer_type>,checker_template_type>::type;
//	using checker_type=Checker<buffer_type>;
    using self_type=session<Protocol,Echo, Consumer, Checker>;
    using error_code_type=boost::system::error_code;
    using pointer_type=std::shared_ptr<self_type>;


	// Create a session from a socket, a consumer (supposing is an external object),  and from a Checker, good that we have a void_t available ....
	session(socket_type socket, consumer_type & consumer, const checker_type & checker=checker_type{}): socket_{std::move(socket)}, checker_{checker},consumer_{consumer}
	{
	}

	~session(){
		std::cout << " ... deleting session ... " << std::endl;
	}

	// Read on start is no parameter
	void start()
	{
		if(started_.load()) return;
		std::cout << " Starting a session ... " << std::endl;
		started_.store(true);
		do_read();
	}

	void stop(){
		if(!started_.load()) return;
        // What's about stop? Do I have to close any service?
		// Shudtown the socket and close it.
		boost::system::error_code ec;
		socket_.shutdown(socket_type::shutdown_both,ec);
		if(ec) std::cout << " Error shutting down socket: " << ec.message() << std::endl;
		socket_.close(ec);
		if(ec) std::cout << " Error closing socket: " << ec.message() << std::endl;
		started_.store(false);
	}

	template<typename Data> void write(const Data * d, std::size_t length)
	{
		write(ba::buffer(d, length));
	}

	template< typename ConstBufferSequence> void write( const ConstBufferSequence & s){
		auto self(this->shared_from_this());
		socket_.async_send(s,
		[this, self](boost::system::error_code ec, std::size_t l /*length*/)
				{
					if (ec)
					{
						std::cout << " Error sending data: " << ec.message() << ", endpoint: " << socket_.local_endpoint()<<", data transfered: "<< l<< std::endl;
					}
				});
	}
	template< typename ConstBufferSequence> void write( const ConstBufferSequence & s, typename socket_type::message_flags flags){
		auto self(this->shared_from_this());
		socket_.async_send(s, flags,
		[this, self](boost::system::error_code ec, std::size_t l /*length*/)
				{
					if (ec)
					{
						std::cout << " Error sending data: " << ec.message() << ", endpoint: " << socket_.local_endpoint()<<", data transfered: "<< l<< std::endl;
					}
				});

	}

	bool is_echo() const {  return Echo; }
	bool is_started() const { return started_.load();}
	socket_type & get_socket(){ return socket_;}

private:
	void do_read()
	{
		// A trick to let the object pointed to, survive to the lambda execution. Otherwise the shared pointer
		// could survive, inside the lambda object, to the object itself
		auto self(this->shared_from_this());
//		socket_.async_read_some(ba::buffer(data_, max_length),
		// If the checker is void, then call the socket member function, otherwise, call the helper function overload with the checker
		if constexpr (std::is_void<checker_template_type>::value) {
			socket_.async_receive(ba::buffer(data_, max_length),
				[this, self](boost::system::error_code ec, std::size_t length)
				{
					if (!ec)
					{
						consumer_(data_,length);
						std::cout << " End of read ... started_:" << std::boolalpha << started_.load()<< std::endl;
						if(started_.load()){
							if constexpr (Echo) {
								do_write(length);
							} else {
								do_read();
							}
						} // if(started_)
					} else {
						std::cout << " An error occurred in the session receive: " << ec.message() << ", endpoint: " << socket_.local_endpoint()<<" \n ... stopping the session \n";
						stop();
					}
				}); // socket_.async_receive
		} else {
			static_assert(!std::is_same<protocol_type,udp>::value," Error in session. UDP socket is not stream oriented. It cannot be used in async_read.");
			ba::async_read(socket_, ba::buffer(data_),
				// The checker functor
    			[this, self] (error_code_type ec, std::size_t length) -> std::size_t
				{
    				if(!ec)
    				{
    					 return read_complete(ec,length);
    				} else {
						std::cout << " An error occurred in the session checker: " << ec.message() << ", endpoint: " << socket_.local_endpoint()<<"\n ... stopping the session \n";
						stop();
    				}
    				return false;
				},
				// The end of read operation handler
				[this, self](boost::system::error_code ec, std::size_t length)
				{
					if (!ec)
					{
						consumer_(data_,length);
						std::cout << " End of read ... started_:" << std::boolalpha << started_.load() << std::endl;
						if(started_.load()){
							if constexpr (Echo) {
								do_write(length);
							} else {
								do_read();
							}
						} // if(started_)
					} else {
						std::cout << " An error occurred in the session end or read : " << ec.message() << ", endpoint: " << socket_.local_endpoint()<<" \n ... stopping the session \n";
						stop();
					}
				}); // ba::asynch_read
		} // if constexpr
	}

	void do_write(std::size_t length)
	{
		auto self(this->shared_from_this());
//		ba::async_write(socket_, ba::buffer(data_, length),
		socket_.async_send(ba::buffer(data_, length),
		[this, self](boost::system::error_code ec, std::size_t /*length*/)
				{
					if (!ec)
					{
						do_read();
					} else {
						std::cout << " An error occurred in the session write: " << ec.message() << " \n ... stopping the session \n";
						stop();
					}
				});
	}

	std::size_t read_complete(const error_code_type & err, std::size_t bytes) {
		// stop reading if an error occurred
		if ( err) return 0;
		// Stop reading if stopped
		std::cout << " Read_complete ... started_:" << std::boolalpha << started_.load() << ", bytes: " << bytes<<std::endl;
		if (!started_.load()) return 0;
		// if no bytes read, then ask for one byte.
		if (bytes == 0) return 1;
		std::size_t bytes_to_read=checker_(data_+(bytes-1));
		std::cout << " Read_complete, bytes_to_read: " << bytes_to_read << std::endl;
		return bytes_to_read;
	}

  socket_type socket_;
  enum { max_length = 1024 };
  buffer_type data_[max_length]={};
  consumer_type & consumer_;
// Use a checker if not void. If a checker is request, then checker_type is used, otherwise a do_nothing class.
//  typename std::conditional<std::is_void<checker_type>::value,do_nothing<buffer_type>,checker_type>::type checker_;
  checker_type checker_;
  std::atomic_bool started_{false};
};


} // namespace networking
} // namespace fcutils


#endif // SESSION_HPP_ 
