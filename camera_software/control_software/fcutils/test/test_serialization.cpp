///
///
/// @brief A new file
///
/// @file test_serializer.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 8 feb 2022 
///

#include <iostream>
#include <string>
#include <cstdint> // for int8_t etc. etc.

#include "string/add_timestamp.hpp"
#include "serialization/serializer.hpp"
#include "serialization/deserializer.hpp"
#include "utility.hpp"

namespace fcser=fcutils::serialization;

struct test{
	// We need to be aligned to the largest type ....
	uint8_t id[2]{0xFAu,0xB0u};
	uint16_t i16[3];
};


int main(){
	fcser::serialize_string a;
	fcser::serialize_repr b;
	fcser::deserialize_string da;
	fcser::deserialize_repr db,ddb;

	std::string val,deser;
	std::cout << " Enter a message: ";
	b.reserve(1024);
	uint16_t count{0u};
	test c;
	while(std::getline(std::cin,val)){
		auto temp=fcutils::string::add_timestamp(val);
		a.clear();
		a.reserve(temp.size());
		a.serialize(temp.c_str(),temp.size());
		deser.reserve(a.data().size());
		da.deserialize(&deser,a.c_str(),a.data().size() );
		std::cout << " Original message: " << val << ", buffer capacity: " << a.data().capacity() << '\n'
				<< " Serialized message: " << a.data() << '\n'
				<< " Deserialized message: " << deser << std::endl;
		b.clear();
		b+="0x";
		b.serialize(temp.c_str(),temp.size());
		std::cout << b.data() << std::endl;
		deser.clear();
		deser.reserve(a.data().size());
		auto pc=deser.data();
		std::size_t j{0};
		// Deserialize using the full string
		db.deserialize(pc,b.c_str(),b.size());
		// Deserialize byte per byte. Remember to skip the first two char from the string size.
		//for(std::size_t i = 2; i<b.data().size(); ++j, i+=2) {
		//	db.deserialize(pc++,b.c_str()+i,2 );
		//}
		std::cout << " Original message: " << val << ", buffer capacity: "<< b.data().capacity() << '\n'
				<< " Serialized message: " << b.data() << ';' << '\n'
				<< " Deserialized message: " << deser.data() << ';'  <<std::endl;
		c.i16[0]=++count;
		c.i16[1]=count*2u;
		c.i16[2]=count*3u;
		b.clear();
		b+="0x";
		b.serialize(&c);
		test d{},dd{};
		auto resdb=db.deserialize(&d,b.c_str(),b.size());
		auto resdd=ddb.deserialize(&dd,b.c_str());
		std::cout << " Original struct: " << std::hex;
		fcutils::print_vals(std::cout, static_cast<uint16_t>(c.id[0]),static_cast<uint16_t>(c.id[1]), c.i16[0], c.i16[1], c.i16[2]);
		std::cout << std::dec <<", buffer capacity: "<< b.data().capacity() << '\n'
				<< " Serialized message: " << b.data() << ';' << '\n'
				<< " Deserialized struct: " << std::hex;
		fcutils::print_vals(std::cout, static_cast<uint16_t>(d.id[0]),static_cast<uint16_t>(d.id[1]), d.i16[0], d.i16[1], d.i16[2]);
		std::cout << " Deserialized struct (size auto deduced): " << std::hex;
			fcutils::print_vals(std::cout, static_cast<uint16_t>(dd.id[0]),static_cast<uint16_t>(dd.id[1]), dd.i16[0],
					dd.i16[1], dd.i16[2]);
		std::cout << resdb  << resdd << std::dec ;
		std::cout << " Enter a message: " << std::endl;
	}
	return 0;
}



