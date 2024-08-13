///
///
/// @brief A new file
///
/// @file serialport.h
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 25 ott 2021 
///

#ifndef SERIAL_PORT_HPP_
#define SERIAL_PORT_HPP_

#include <boost/asio/serial_port.hpp>
#include <boost/asio/io_service.hpp>
#include <serialport/serialopts.hpp>
#include <string>

// The serial port options

namespace fcutils{
namespace serialport{

namespace ba=boost::asio;

struct serial_port{

	using options_type=fcutils::serialport::options;

	serial_port(ba::io_service &s, const std::string &n): name(n), service(s), port(s,n) {
	  fcutils::serialport::get_port_opts(port,origoptions);
	};

	serial_port(ba::io_service &s, const std::string &n, const options_type &opts): name(n), service(s), port(s,n) {
		  fcutils::serialport::get_port_opts(port,origoptions);
		  set_options(opts);
	};
	// Copy constructor & assignment must be deleted because ASIO serial port cannot be copied
	serial_port(const serial_port & s) = delete;
	serial_port & operator=(const serial_port & s) = delete;
	// Move constructor & assignment (No need to move the service - is a reference- and the origoptions - no deep copies involved-)
	serial_port(serial_port && s) noexcept : port{std::move(s.port)}, name{std::move(s.name)}, service{s.service}, origoptions{s.origoptions}{};
	serial_port & operator=(serial_port && s) =delete;
	// Do I need to close the port or cancel operation?
	~serial_port(){ set_options(origoptions);}

  	template <class O> void get_option(O & opt) const { port.get_option(opt); }
  	template <class O> void set_option(O & opt){ port.set_option(opt); }
  	options_type get_options() const { return fcutils::serialport::get_port_opts(port);}
  	void set_options(const options_type & opts ) {fcutils::serialport::set_port_opts(port,opts); }
// Data members
  	ba::serial_port port;
	std::string name;
	ba::io_service &service;
	options_type origoptions;
};

std::ostream & operator <<(std::ostream & o, const serial_port &sp){
	serial_port::options_type temp{sp.get_options()};
	o << "Serial Port on: " << sp.name << ", actual options. \n"
			<<	temp;
	return o;
}

} // namespace serialport
} // namsepace fcutils

#endif // SERIAL_PORT_HPP_
