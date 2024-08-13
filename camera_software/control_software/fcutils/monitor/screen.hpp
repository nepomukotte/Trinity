///
///
/// @brief A new file
///
/// @file screen.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 26 mar 2022 
///

#ifndef SCREEN_HPP_
#define SCREEN_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <optional>
#include <cmath>
#include <tuple>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "fcutils/terminal/ansi.hpp"
#include "fcutils/monitor/favorite.hpp"
#include "fcutils/monitor/xmlconst.hpp"
#include "fcutils/monitor/color_map.hpp"
#include "fcutils/monitor/monitor_const.hpp"


namespace fcutils{
namespace monitor{

template <typename Value, typename Data> struct screen{
	using favorite_type 		= favorite<Value,Data>;
	using favorite_vector_type	= std::vector<favorite_type>;
	using screensize_type		= std::pair<size_t,size_t>;
//	using tree_type 			= favorite_vector_type::value_type::tree_type;
	using tree_type 			= boost::property_tree::ptree;
	using data_type				= typename favorite_vector_type::value_type::data_type;
	using favorite_name_type	= std::string;
	using favorite_filename_type= std::string;
	using favorite_list_type	= std::vector< std::pair<favorite_name_type, favorite_filename_type> >;
	using color_modifier_type	= fcutils::terminal::ansi::color;

	std::string title{"None"};
	std::string info{"No Info"};
	size_t columns{0};
	size_t lines{0};
    color_modifier_type textcolor{fcutils::terminal::ansi::color::FG_DEFAULT};
    color_modifier_type bgcolor{fcutils::terminal::ansi::color::BG_GREEN};

	favorite_vector_type favorites;
	screensize_type screensize;
    favorite_list_type favlistnames;

    void print(std::ostream & os=std::cout) const {
		os << " Monitor name: " << title <<", info: "<< info <<'\n'
		   << " Columns: " << columns <<", lines:" << lines << '\n'
		   << " Favorites names and filenames found: "<< favlistnames.size() << ", lists: \n";
		for( const auto & [fav_name, fav_file]: favlistnames){
			os << " Name: " << fav_name <<", filename: " << fav_file << '\n';
		}
		os << std::endl;
    }

    void print_all(std::ostream & os=std::cout) const{
    	print(os);
    	for( const auto & v: favorites){
    		v.print_all(os);
    	}
    }

    bool fetch(const std::string & filename){
    	tree_type pt;
    	read_xml(filename, pt);
    	return fetch(pt);
    }

    bool fetch( tree_type & pt){
    	bool success=true;
    	auto ptmonitor=pt.get_child_optional(xml::monitor::MONITORROOT);
    	if(ptmonitor) {
    		title=ptmonitor->get_child(xml::monitor::field::TITLEOPT).get_value("notitle");
    		info=ptmonitor->get_child(xml::monitor::field::INFOOPT).get_value("noinfo");
    		columns=ptmonitor->get_child(xml::monitor::field::COUMNOPT).get_value(XTERMCOLUMN);
    		auto ptvar=ptmonitor->get_child_optional(xml::slowcontrol::field::TEXTCOLOROPT);
    		if(ptvar) {
    			std::string temp=ptvar->get_value(xml::slowcontrol::favorite::value::SRED);
    			textcolor=string_to_object(temp, STRING_TO_COLOR);
    		}
    		ptvar.reset();
    		ptvar=ptmonitor->get_child_optional(xml::slowcontrol::field::BACKGROUNDCOLOROPT);
    		if(ptvar) {
    			std::string temp=ptvar->get_value(xml::slowcontrol::favorite::value::SRED);
    			bgcolor=string_to_object(temp, STRING_TO_BG_COLOR);
    		}

    		auto ptfavorites=pt.get_child_optional(xml::monitor::FAVORITESROOT);
        	if(ptfavorites){
        		for( const auto & [name,ptree]: ptfavorites.value() ){
        			favlistnames.emplace_back(
        					ptree.get_child(xml::monitor::field::FAVORITENAMEOPT).get_value("noname"),
							ptree.get_child(xml::monitor::field::FAVORITEFILEOPT).get_value("nofile")
        			);
        		}
    		} // if(ptfavorites)
    	} // if(ptmonitor)
		return success;
	}

    bool fetch_all(const std::string & filename){
    	tree_type pt;
    	read_xml(filename, pt);
    	return fetch_all(pt);
    }

	bool fetch_all( tree_type & pt, const std::string &favname=xml::slowcontrol::favorite::value::DEFAULTNAME){
		bool success=false;
		success=fetch(pt);
		success=fetch_favorites();
		return success;
	}

	bool fetch_favorite( const favorite_name_type & favname, const favorite_filename_type & filename ){
		bool success=true;
		tree_type pt;
		read_xml(filename,pt);
		// Here we need the root property tree because we need to fetch also the variables tree in the favorite::fetch_all
		favorite_type temp;
		success=temp.fetch_all(pt,favname);
		if(success) favorites.push_back(temp);
		return success;

	}

	bool fetch_favorites(){
		bool success=true,temp;
		for( const auto & [name,file]: favlistnames){
			temp=fetch_favorite(name,file);
			success=(temp)?success:temp;
		}
		return success;
	}

	screensize_type screen_size(size_t extra_lines=0) {
		// calculate maximum label size
		screensize_type::first_type maxlabel=0;
		// Calculate the lines number. One line for the monitor title, one for each favorite title, 2 for each value and 1 for the last messages
		// Add to chars to the maximumlabel size.
		// Keep two extra lines for some printout
		lines=2;
		// Add the extra_lines defined by the user
		lines+=extra_lines;
		// To calculate the size compare the label width and the field width
		for( const auto & v : favorites){
			maxlabel=(maxlabel<v.width)?v.width:maxlabel;
			for ( const auto & vv : v.values ){
				maxlabel=(maxlabel<vv.label.size())?vv.label.size():maxlabel;
			}
			// Calculate the lines number. One line for the monitor title, one for each favorite title, 2 for each value and 1 for the last messages
			lines+=2*std::ceil(static_cast<double>(v.values.size())/static_cast<double>(columns))+1;
		}
		screensize=screensize_type( (maxlabel+2)*columns,lines);
		return screensize;
	}


    void set_values( data_type *d) {
    	for( auto & v: favorites){
    		v.set_values(d);
    	}
    }

    void init( data_type d=data_type{0} ){
    	for( auto & fav: favorites){
    		fav=d;
    	}
    }

}; // screen

// Helper functions to print parts of the screen collection

template <typename Screen> void list_values( const Screen & scmon, std::ostream &os=std::cout){
	for( auto & fav: scmon.favorites){
		bool ifirst=true;
		char separator=' ';
		for( auto & vv: fav.values){
			separator=ifirst?' ':',';
			ifirst=false;
			os << separator << vv ;
		}
	}
	os << std::endl;

}

template < typename Screen> void list_labels( const Screen & scmon, std::ostream &os=std::cout){
	std::ios_base::fmtflags ff;
	// Backup the ostream format
	ff = os.flags();
	for(auto & fav: scmon.favorites) {
		bool ifirst=true;
		char separator=' ';
		for( auto & vv: fav.values){
			os.setf (vv.outfmt , std::ios::basefield);
			separator=ifirst?' ':',';
			ifirst=false;
			os << separator << vv.label ;
		}
	}
	os.flags(ff);
	os << std::endl;
}

template < typename Screen> void adapt_terminal( Screen & scmon, size_t extra_lines=0, std::ostream & os=std::cout){
	// Calculate the screen size
	auto [columns,lines]= scmon.screen_size(extra_lines);
	// Reset the screen
	fcutils::terminal::ansi::resize_clear(lines,columns,os);
	// Position the cursor at the begining
	fcutils::terminal::ansi::gotoxy(0ul,0ul);
}

template <typename ValueCollection>
void print_values( const ValueCollection & values, size_t column, size_t precision,
		size_t width, std::ostream & os=std::cout){
	int stop=0,start=stop;
	std::stringstream fmtlabel;
	fmtlabel << "%=" << width <<'s';
    auto m=fcutils::terminal::ansi::color::BG_DEFAULT;
    boost::format fmtvalue;
	for(auto & v: values){
		++stop;
		std::cout << fcutils::terminal::ansi::format::TXT_BOLD << boost::format(fmtlabel.str().c_str()) % v.label ;
		if( ( (stop % column)==0 ) || stop==values.size() ) {
			std::cout << fcutils::terminal::ansi::format::TXT_RESET << '\n';
			for(int i=start; i< stop; ++i){
				// Modify the BOOST format object according to value output format
				std::cout << values[i].to_format(fmtvalue,width,precision);
			}
		std::cout << std::endl;
		start=stop;
		}
	}
}

template < typename Favorite>
void print_favorites( const Favorite &fav, size_t column, std::size_t xscreensize, //=::monitor::XTERMCOLUMN*::monitor::FIELDWIDTH ,
		std::ostream & os=std::cout){
	std::stringstream fmt;
	fmt << "%=" << xscreensize <<'s';
    std::cout << fav.bgcolor << fcutils::terminal::ansi::format::TXT_BOLD << fav.textcolor << boost::format(fmt.str().c_str()) % fav.name
    		<< fcutils::terminal::ansi::format::TXT_RESET << std::endl;
    print_values(fav.values, column, fav.precision, static_cast<size_t>(xscreensize/column),os);
}

template <typename Value, typename Data> std::ostream & operator<<(std::ostream & os, const screen<Value,Data> & scmon ){
	std::stringstream fmt;
	fmt << "%=" << (scmon.screensize.first) <<'s';
	std::cout << fcutils::terminal::ansi::format::TXT_RESET << std::endl;
    std::cout << scmon.bgcolor << fcutils::terminal::ansi::format::TXT_BOLD << scmon.textcolor << boost::format(fmt.str().c_str()) % scmon.title
    		<< fcutils::terminal::ansi::format::TXT_RESET << std::endl;
    for(auto & fav: scmon.favorites) {
		print_favorites(fav,scmon.columns, scmon.screensize.first, os);
	}
    return os;
}

} // namespace monitor
} // namespace fcutils




#endif // SCREEN_HPP_ 
