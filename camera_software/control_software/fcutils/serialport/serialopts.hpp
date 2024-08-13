///
///
/// @brief A new file
///
/// @file serialopts.h
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 25 ott 2021 
///

#ifndef SERIALOPTS_HPP_
#define SERIALOPTS_HPP_

#include <iostream>
#include <string>
// Header for asio serial port
#include <boost/asio/serial_port.hpp>
#include <boost/asio/serial_port_base.hpp>

namespace fcutils{
namespace serialport{

namespace ba=boost::asio;

// Serial options to be used to set up the serial port

struct options{
// Just rely on the default constructors
	ba::serial_port_base::baud_rate      rate{};
	ba::serial_port_base::flow_control   control{};
	ba::serial_port_base::parity         parity{};
	ba::serial_port_base::stop_bits      stop{};
	ba::serial_port_base::character_size size{};
};

// setters and getters
void set_opts(options & s, ba::serial_port_base::baud_rate    &b){ s.rate=b; }
void set_opts(options & s, ba::serial_port_base::flow_control &f){ s.control=f; }
void set_opts(options & s, ba::serial_port_base::parity       &p){ s.parity=p; }
void set_opts(options & s, ba::serial_port_base::stop_bits    &sb){ s.stop = sb;}
void set_opts(options & s, ba::serial_port_base::character_size &c){ s.size=c;}

options get_port_opts(const boost::asio::serial_port & p){
	options temp;
	p.get_option(temp.rate);
	p.get_option(temp.control);
	p.get_option(temp.stop);
	p.get_option(temp.parity);
	p.get_option(temp.size);
	return temp;
}

void get_port_opts(const boost::asio::serial_port &p, options &s){
	s=get_port_opts(p);
}

void set_port_opts(boost::asio::serial_port &p, const options &s){
	p.set_option(s.rate);
	p.set_option(s.control);
	p.set_option(s.stop);
	p.set_option(s.parity);
	p.set_option(s.size);
}


// Few helper function to convert the enum into strings. It is just for printing purpose
const std::string value_name(const ba::serial_port_base::flow_control & p){
	switch(p.value()){
	case ba::serial_port_base::flow_control::none  :
		return "none";
	case ba::serial_port_base::flow_control::software  :
		return "software";
	case ba::serial_port_base::flow_control::hardware  :
		return "hardware";
	}
	return "ERROR";
}

const std::string value_name(const ba::serial_port_base::parity & p){
	switch(p.value()){
	case ba::serial_port_base::parity::none  :
		return "none";
	case ba::serial_port_base::parity::odd  :
		return "odd";
	case ba::serial_port_base::parity::even  :
		return "even";
	}
	return "ERROR";
}

const std::string value_name(const ba::serial_port_base::stop_bits & p){
	switch(p.value()){
	case ba::serial_port_base::stop_bits::one  :
		return "one";
	case ba::serial_port_base::stop_bits::onepointfive  :
		return "onepointfive";
	case ba::serial_port_base::stop_bits::two  :
		return "two";
	}
	return "ERROR";
}

ba::serial_port_base::flow_control string2serialopt( const std::string & p, ba::serial_port_base::flow_control & d ){
	ba::serial_port_base::flow_control::type temp;
	if(p== "none") temp=ba::serial_port_base::flow_control::none;
	if(p== "software") temp=ba::serial_port_base::flow_control::software;
	if(p== "hardware") temp=ba::serial_port_base::flow_control::hardware;
	return ba::serial_port_base::flow_control(temp);
}

ba::serial_port_base::parity string2serialopt( const std::string & p, ba::serial_port_base::parity & d){
	ba::serial_port_base::parity::type temp;
	if(p=="none") temp= ba::serial_port_base::parity::none;
	if(p=="odd" ) temp= ba::serial_port_base::parity::odd;
	if(p=="even") temp= ba::serial_port_base::parity::even;
	return ba::serial_port_base::parity(temp);
}

ba::serial_port_base::stop_bits string2serialopt( const std::string & p, ba::serial_port_base::stop_bits & d){
	ba::serial_port_base::stop_bits::type temp;
	if(p=="one") 		temp=ba::serial_port_base::stop_bits::one;
	if(p=="onpointfive")temp=ba::serial_port_base::stop_bits::onepointfive;
	if(p=="two") temp=ba::serial_port_base::stop_bits::two;
	return ba::serial_port_base::stop_bits(temp);
}

ba::serial_port_base::baud_rate string2serialopt (const std::string &p, ba::serial_port_base::baud_rate & d ){
	int v=std::stoi(p);
	return ba::serial_port_base::baud_rate(v);
}

ba::serial_port_base::character_size string2serialopt (const std::string &p, ba::serial_port_base::character_size & d ){
	int v=std::stoi(p);
	return ba::serial_port_base::character_size(v);
}

// A stream operator for the serial options
std::ostream & operator<<(std::ostream &o, const options & s){
	o<< "Baud Rate     : " << s.rate.value() << "\n"
	  <<"Flow Control  : " << value_name(s.control) << "\n"
	  <<"Parity        : " << value_name(s.parity) << "\n"
	  <<"Stop Bits     : " << value_name(s.stop) << "\n"
	  <<"Character Size: " << s.size.value() << "\n";
	return o;
}

} // namespace serialport
} // namespace fcutils


#endif // SERIALOPTS_HPP_ 
