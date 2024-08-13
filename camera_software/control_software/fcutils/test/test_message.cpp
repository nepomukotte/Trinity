///
///
/// @brief A new file
///
/// @file test_message.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 9 mar 2022 
///
#include <iostream>
#include <cstddef>
#include <chrono>
#include <cstring>
#include <type_traits>
#include "fcutils/utility/utility.hpp"
#include "fcutils/message/time_ht.hpp"
#include "fcutils/message/general_ht.hpp"
#include "fcutils/message/crc_ht.hpp"

template <typename Data> struct base_consumer{
	using data_type=Data;
	using data_ptr_type=Data *;
	virtual void consume(const data_type *, std::size_t l)=0;
	virtual void print_stats(std::ostream & o, const std::string &prefix=""){
		o << prefix << " Nothing to be printed \n";
	}
	virtual void reset(){
		std::cout << "I am resetting ... " << std::endl;
	}
};

template< typename Data> struct print: public base_consumer<Data>{
	void consume(const  typename base_consumer<Data>::data_type * d, std::size_t l){
		operator()(d,l);
	}
	void operator()(const typename base_consumer<Data>::data_type * d, std::size_t l){
		std::cout << " Bytes received: " << l <<", values:\n";
		std::ostreambuf_iterator<typename base_consumer<Data>::data_type> oo (std::cout.rdbuf()); // output iterator for cout
		std::copy(d,d+l,oo);
		std::cout << std::endl;
	}
};

template <typename Data> struct print_byte: public base_consumer<Data>{
	void consume(const typename base_consumer<Data>::data_type * d, std::size_t l){
		operator()(d,l);
	}
	void operator()(const typename base_consumer<Data>::data_type * d, std::size_t l){
		std::cout << " Bytes received: " << l <<", values:\n";
		for(int i=0; i<l; ++i) std::cout << std::hex<< std::to_integer<short>(static_cast<std::byte>(static_cast<unsigned char>(*(d+i)))) <<' ';
		std::cout << std::dec << std::endl;
	}
};

struct count_byte{
	using Data = unsigned char;
	unsigned short consume(const Data * d, std::size_t l){
		std::cout << " Consuming ..." << std::endl;
		return operator()(d,l);
	}
	unsigned short operator()(const Data * d, std::size_t l){
		std::size_t s=sizeof(Data);
		fcutils::print_vals(std::cout,s,l,s*l);
		return static_cast<unsigned short>(s*l);
	}
};

int main(){
	using header_type=fcutils::message::time_ht<unsigned char, 15,fcutils::message::TIME_HT_OPTION::AT_END>;
	header_type a={0xab,0xcd,0xef};
//	std::array<unsigned char, 5> a ={0xab,0xcd,0xef};
	std::cout << std::boolalpha << " Header callable: " << std::is_invocable_v<
			decltype(static_cast< void (header_type::*)(const header_type::data_type *, std::size_t) >(&header_type::consume)), header_type *, const header_type::data_type *, std::size_t
			> << std::hex << std::endl;
	std::cout << std::boolalpha << " Header callable: " << std::is_invocable_v<
			decltype(static_cast< void (header_type::*)(const header_type::data_type *, std::size_t) >(&header_type::operator())), header_type *, const header_type::data_type *, std::size_t
			> << std::hex << std::endl;
	print_byte<unsigned char> c;
	a();
    fcutils::print_vals<':'>(std::cout," Offset ",a.time_offset(),sizeof(std::time_t));
	c.consume(a.data(),a.size());
	std::cout << std::ctime(reinterpret_cast<std::time_t*>(a.data()+a.time_offset())) << std::endl;
	std::time_t t=a;
	std::cout << std::ctime(&t) << std::endl;
	fcutils::message::general_ht<unsigned char, 16, count_byte, fcutils::message::GENERAL_HT_OPTION::AT_END> b={0xab,0xcd,0xef};
	b();
	c.consume(b.data(),b.size());
    fcutils::print_vals<':'>(std::cout," Offset ",b.value_offset(),sizeof(unsigned short));

    fcutils::message::crc16_ht<unsigned char> d{0xFA,0xB0};
    d(a.data(),a.size());
    c.consume(d.data(),d.size());
    fcutils::print_vals<':'>(std::cout," Offset ",d.value_offset(),sizeof(unsigned short)," Latest ", d.latest_value(),*reinterpret_cast<uint16_t*>(d.data()+d.value_offset()));
//	std::cout << std::dec << std::endl;

	return 0;

}


