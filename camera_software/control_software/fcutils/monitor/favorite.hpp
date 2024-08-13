///
///
/// @brief A new file
///
/// @file favorite.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 26 mar 2022 
///

#ifndef FAVORITE_HPP_
#define FAVORITE_HPP_
#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <algorithm>

#include <boost/format.hpp>
#include <boost/property_tree/ptree.hpp>

#include "fcutils/terminal/ansi.hpp"
#include "fcutils/monitor/monitor_const.hpp"
#include "fcutils/monitor/xmlconst.hpp"
#include "fcutils/monitor/color_map.hpp"
#include "fcutils/property_tree/property_tree.hpp"

namespace fcutils{
namespace monitor{

template <typename Value, typename Data = unsigned char> struct favorite{
	using self_type				= favorite;
	using value_type 			= Value;
	using tree_type				= typename value_type::tree_type;
//	using data_type				= typename value_type::data_type;
	using data_type				= Data;
	using value_vector_type		= std::vector<value_type>;
	using favorite_vector_type	= std::vector< std::string> ;
	using format_type			= boost::format;
	using color_type			= fcutils::terminal::ansi::color;

	std::string name{"Noname"};
	std::string info{"No info"};
	size_t precision{};
	size_t width{};
	value_vector_type values;
    favorite_vector_type favorites;
    format_type fmt;
    color_type textcolor{fcutils::terminal::ansi::color::FG_BLUE};
    color_type bgcolor{fcutils::terminal::ansi::color::BG_GREEN};

	template<typename Val> self_type & operator=(const Val val){
    	for( auto & v: values){
    		v=val;
    	}
	}

    void print(std::ostream & os=std::cout) const {
		os << " Favorite list name: " << name <<", info: "<< info <<'\n';
		os << " Floating point precision: " << precision <<", field width :" << width << '\n'; // ", boost format object: " << fmt.str() <<'\n';
		os << " Favorite variables: "  << favorites.size() << ", names: \n";
		std::copy(favorites.begin(),
				  favorites.end(),
				  std::ostream_iterator<favorite_vector_type::value_type>(os," "));
		os << std::endl;
    }

    void print_all(std::ostream & os=std::cout) const{
    	print(os);
    	for( const auto & v: values){
    		v.print(os);
    	}
    }

    void set_values(data_type * d) {
    	for( auto & v: values){
    		v.set_value(d);
    	}
    }

    bool fetch( tree_type & pt, const std::string & favname=xml::slowcontrol::favorite::value::DEFAULTNAME){
    	bool success=true;
    	auto ptfavorites=fcutils::property_tree::find_by_value_optional(pt,favname);
    	if(ptfavorites) {
    		name=ptfavorites->get_child(xml::slowcontrol::favorite::field::NAMEOPT).get_value("noname");
    		info=ptfavorites->get_child(xml::slowcontrol::favorite::field::INFOOPT).get_value("noinfo");
    		precision=ptfavorites->get_child(xml::slowcontrol::favorite::field::PRECISIONOPT).get_value(PRECISION);
    		width=ptfavorites->get_child(xml::slowcontrol::favorite::field::WIDTHOPT).get_value(FIELDWIDTH);
    		auto ptvar=ptfavorites->get_child_optional(xml::slowcontrol::favorite::field::TEXTCOLOROPT);
    		if(ptvar) {
    			std::string temp=ptvar->get_value(xml::slowcontrol::favorite::value::SRED);
    			textcolor=string_to_object(temp, STRING_TO_COLOR);
    		}
    		ptvar.reset();
    		ptvar=ptfavorites->get_child_optional(xml::slowcontrol::favorite::field::BACKGROUNDCOLOROPT);
    		if(ptvar) {
    			std::string temp=ptvar->get_value(xml::slowcontrol::favorite::value::SRED);
    			bgcolor=string_to_object(temp, STRING_TO_BG_COLOR);
    		}
    		std::stringstream ss;
    		ss<<'%' << width << '.' << precision << 'f';
    		fmt.clear_binds();
    		fmt=boost::format(ss.str().c_str());
    		for( auto & [key,fav]: *ptfavorites ){
    			if(auto nodename=fav.get_child_optional(xml::slowcontrol::favorite::field::NAMEOPT)) {
    							// favorite.name found
    				favorites.push_back(nodename->get_value("no name"));
    			}

    		}
    	} // if(ptfavorites)
		return success;
	}


	bool fetch_values( tree_type & pt, const std::string & valroot=xml::slowcontrol::VARIABLESROOT){
		bool success=true;
		auto ptvariables=pt.get_child_optional(valroot);
		if(ptvariables) {
			for( auto & v: favorites){
				boost::optional<tree_type& > ptvalue=fcutils::property_tree::find_by_value_optional(*ptvariables,v);
				if(ptvalue){
					values.emplace_back().fetch(*ptvalue);
				}
			}
		}// if(ptvariables)
		return success;
	}

	bool fetch_all( tree_type & pt, const std::string &favname=xml::slowcontrol::favorite::value::DEFAULTNAME){
		bool success=false;
		success=fetch(pt,favname);
		success=success&&fetch_values(pt);
		return success;
	}
};


} // namespace monitor
} // namespace fcutils



#endif // FAVORITE_HPP_ 
