///
///
/// @brief A new file
///
/// @file find_by_value_optional.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 26 mar 2022 
///

#ifndef FIND_BY_VALUE_OPTIONAL_HPP_
#define FIND_BY_VALUE_OPTIONAL_HPP_
#include <iostream>

#include <boost/property_tree/ptree.hpp>

namespace {
namespace bpt=boost::property_tree;

}
namespace fcutils{
namespace property_tree{

// Search using a key for a given data value
boost::optional<bpt::ptree &> find_by_value_optional( bpt::ptree &pt, const bpt::ptree::data_type & data ){
	boost::optional<bpt::ptree &> temp;
	bpt::ptree::data_type tempdata;
	for( auto & [key,subnode]: pt) {
		if(!subnode.empty()){
			// The node has a subnode launch the search on the subnode
			// If has been found, then the optional assignment return true
			temp=find_by_value_optional(subnode, data);
			// If a valid tree has been returned, then return
			if(temp) return temp;
		}
		// Here the node ha no subnode, just compare the node value
		auto tempdata=subnode.data();
		if(tempdata==data){
			std::cout << " Value: " << data << ", found in node: " << key << std::endl;
			temp = pt;
			return temp;
//			return pt;
		}
	} // BOOST_FOREACH
	return temp;
}
//		}
//	}
//	BOOST_FOREACH( bpt::ptree::value_type & v, pt){
//		if(!v.second.empty()){
//			// The node has a subnode launch the search on the subnode
//			// If has been found, then the optional assignment return true
//			temp=find_by_value_optional(v.second, data);
//			// If a valid tree has been returned, then return
//			if(temp) return temp;
//		}
//		// Here the node ha no subnode, just compare the node value
//		tempdata=v.second.data();
//		if(tempdata==data){
//			std::cout << " Value: " << data << ", found in node: " << v.first.data() << std::endl;
//			return pt;
//		}
//	} // BOOST_FOREACH
//	return temp;
//}

} // namespace property_tree
} // fcutils



#endif // FIND_BY_VALUE_OPTIONAL_HPP_ 
