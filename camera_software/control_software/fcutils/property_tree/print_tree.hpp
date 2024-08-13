///
///
/// @brief A new file
///
/// @file print_tree.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 26 mar 2022 
///

#ifndef PRINT_TREE_HPP_
#define PRINT_TREE_HPP_
#include <iostream>

#include <boost/property_tree/ptree.hpp>

namespace {
namespace bpt=boost::property_tree;
}
namespace fcutils{
namespace ptree{

void printtree(bpt::ptree &pt){
	for( auto & [key,tree]: pt){
		std::cout << "node found: " << key.data();
		if(!tree.empty()) {
			std::cout << " subvalue found for node: " << key.data() << std::endl;
			printtree(tree);
		} else {
			std::cout << ", Value: " << tree.data() << std::endl;
		}

	}
//	BOOST_FOREACH( bpt::ptree::value_type & v, pt){
//		std::cout << "node found: " << v.first.data();
//		if(!v.second.empty()) {
//			std::cout << " subvalue found for node: " << v.first.data() << std::endl;
//			printtree(v.second);
//		} else {
//			std::cout << ", Value: " << v.second.data() << std::endl;
//		}
//
//	}
}

} // namespace ptree


} // namespace fcutils




#endif // PRINT_TREE_HPP_ 
