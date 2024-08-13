///
///
/// @brief A new file
///
/// @file talk_to_serial.h
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 26 ott 2021 
///

#ifndef TALK_TO_SERIAL_HPP_
#define TALK_TO_SERIAL_HPP_

#include <boost/asio.hpp>
#include <boost/core/noncopyable.hpp>
#include <serialport/serial_port.hpp>
#include <string>
#include <memory> // std::enable_shared_from_this, std::shared_ptr
#include <algorithm> // std::find
#include <utility> // std::forward


namespace fcutils{
namespace serialport{

namespace ba=boost::asio;

// Basic consumer
template<class Data> struct do_nothing{
    void operator()(Data *, std::size_t){};
};

template<class Data> struct print{
    void operator()(Data* d, std::size_t l){
        if(l>0){
            //std::cout << " Bytes received: " << l <<", values:\n";
            std::ostreambuf_iterator<Data> oo (std::cout.rdbuf()); // output iterator for cout
            std::copy(d,d+l,oo);
            std::cout << std::endl;
        }
    }
};
// Basic checker
// Here we suppose Data is char ....
template<class Data> struct check_eol{
    size_t operator()(Data* d, size_t bytes) {
        // We read one byte at the time until we get '\n'
        // an alternative way checking the full received array
        //      bool found = ( std::find(d, d + bytes, '\n') < (d + bytes) );
        bool found=(*(d+bytes-1)=='\n');
        return found ? 0 : 1;
    }
};


template<class Data=char, template<typename> class Consumer=print, template<typename> class Checker=check_eol >
class talk_to_serial : public std::enable_shared_from_this< talk_to_serial<Data,Consumer,Checker> >, boost::noncopyable {
public:
    using data_type         =Data;
    using consumer_type     =Consumer<data_type>;
    using checker_type      =Checker<data_type>;
    using self_type         =talk_to_serial;
    using error_code_type   =boost::system::error_code;
    using pointer_type      =std::shared_ptr<self_type>;
    using serial_port_type  =fcutils::serialport::serial_port;

    void start_read() {
        started_ = true;
        do_read();
    }

//    static pointer_type get(ba::io_service &s, const std::string &pn) {
//        pointer_type temp(new talk_to_serial(s,pn) );
//        return temp;
//    }

    // Use universal references to rule them all ... Remember to always forward universal references
    template<typename ... Arg >   static pointer_type get(Arg && ... args) {
        pointer_type temp(new talk_to_serial(std::forward<Arg>(args)...) );
        return temp;
    }

    void stop_read() {
        started_ = false;
        // What else? Do I have to close anything?
    }

    serial_port_type  & get_port() { return serialport_;}

    void write(const data_type * data, std::size_t length) {
        // Send a byte at the time to avoid data pileup
        for (int i=0; i<length;++i) ba::write(serialport_.port,ba::buffer(data+i,1));
    }

private:
    talk_to_serial(ba::io_service &s, const std::string &pn, const consumer_type & con, const checker_type & che) :
        serialport_(s, pn), started_(false), consumer_{con}, checker_{che} {};
    talk_to_serial(ba::io_service &s, const std::string &pn,const consumer_type & con) :
            serialport_(s, pn), started_(false), consumer_{con} {};
    talk_to_serial(ba::io_service &s, const std::string &pn, const checker_type & che) :
                serialport_(s, pn), started_(false), checker_{che} {};
    talk_to_serial(ba::io_service &s, const std::string &pn) :
            serialport_(s, pn), started_(false) {}

    void on_read(const error_code_type & err, size_t bytes) {
        if(err){
            std::cout << " Error reading serial port " << serialport_.name << ", bytes red: " << bytes << ", error message: \n"
                    << err.message() << std::endl;
        } else {
            // Do something with the data ...
            consumer_(read_buffer_, bytes);
        } // if (err)
        // Read another message if not stopped
        if(started_) do_read();
    }

    void do_read() {
        // Asynch read a message
        auto self(this->shared_from_this());
        ba::async_read(serialport_.port, ba::buffer(read_buffer_),
                // The Checker implement the protocol checking and the amount of bytes to read further
                [this, self] (error_code_type ec, std::size_t length) -> std::size_t
                {
                    if(!ec)
                    {
                         return read_complete(ec,length);
                    }
                    return false;
                },
                // The Handler called at the end of the readout
                [this,self] (error_code_type ec, std::size_t length)
                {
                    on_read(ec,length);
                });

    }

    std::size_t read_complete(const error_code_type & err, std::size_t bytes) {
        if ( err) return 0;
        // Stop reading if stopped
        if (!started_) return 0;
// if no bytes read, then ask for one byte.
        if (bytes == 0) return 1;
        std::size_t bytes_to_read=checker_(read_buffer_,bytes);
        return bytes_to_read;
     }

private:
    serial_port_type serialport_;
    consumer_type consumer_;
    checker_type checker_;
    enum { max_msg = 1024 };
    data_type read_buffer_[max_msg];
    bool started_;
};

} // namespace serialport
} // namespace fcutils


#endif // TALK_TO_SERIAL_HPP_
