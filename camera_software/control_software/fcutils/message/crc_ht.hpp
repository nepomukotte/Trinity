///
///
/// @brief An header or trailer with a CRC. Based on the 16 bit CRC used for SPB1
///
/// @file crc_ht.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 13 mar 2022 
///

#ifndef CRC_HT_HPP_
#define CRC_HT_HPP_

#include <cstdint> // for intX_t, uintX_t, and fast and least version
#include <algorithm> // for std::for_each
#include <boost/crc.hpp>

#include "fcutils/message/general_ht.hpp"


namespace fcutils{
namespace message{

// The CRC related const
namespace crc {
	typedef unsigned char byte;
	const unsigned int POLY=0x5935;
	const unsigned int START=0x0;
	const unsigned int FINAL=0x0;
	const int SIZE=16;
	const bool REFLECTIN=false;
	const bool REFLECTOUT=false;
    using  boost_crc16_type=boost::crc_optimal<SIZE, POLY, START,
    	FINAL, REFLECTIN, REFLECTOUT>;
    template <typename Data> struct crc16{
    	using data_type=Data;
    	boost_crc16_type crc;
    	uint16_t operator()(const data_type * c, std::size_t length){
    		crc.reset();
    		return std::for_each( c, c + length , crc )();

    	}
    	uint16_t consume(const data_type * c, std::size_t length){
    		return operator()(c,length);
    	}
    };
} // namespace CRC

template<typename Data> using crc16_ht = general_ht<Data, 4, crc::crc16<Data>, GENERAL_HT_OPTION::AT_END>;

} // message
} // fcutils




#endif // CRC_HT_HPP_ 
