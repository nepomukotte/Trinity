///
///
/// @brief A new file
///
/// @file test_monitor.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 29 mar 2022 
///
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "fcutils/monitor/value.hpp"
#include "fcutils/monitor/screen.hpp"

namespace po=boost::program_options;
namespace pt=boost::property_tree;

// Search using a key for a given data value
boost::optional<pt::ptree &> find_by_value_optional( pt::ptree &ptree, const pt::ptree::data_type & data ){
	boost::optional<pt::ptree &> temp;
	pt::ptree::data_type tempdata;
//	BOOST_FOREACH( bpt::ptree::value_type & v, ptree){
	for( auto & v: ptree) {
		if(!v.second.empty()){
			// The node has a subnode launch the search on the subnode
			// If has been found, then the optional assignment return true
			temp=find_by_value_optional(v.second, data);
			// If a valid tree has been returned, then return
			if(temp) return temp;
		}
		// Here the node ha no subnode, just compare the node value
		tempdata=v.second.data();
		if(tempdata==data){
			std::cout << " Value: " << data << ", found in node: " << v.first.data() << std::endl;
			return ptree;
		}
	} // BOOST_FOREACH
	return temp;
}

pt::ptree find_by_value( pt::ptree &ptree, const pt::ptree::data_type & data ){
	auto temp=find_by_value_optional( ptree, data );
	if(temp) return *temp;
	return pt::ptree{};
}

void printtree(pt::ptree &ptree){
//	BOOST_FOREACH( pt::ptree::value_type & v, ptree){
	for ( auto & v: ptree){
		std::cout << "node found: " << v.first.data();
		if(!v.second.empty()) {
			std::cout << " subvalue found for node: " << v.first.data() << std::endl;
			printtree(v.second);
		} else {
			std::cout << ", Value: " << v.second.data() << std::endl;
		}

	}
}

int main(int ac, char* av[])
{
	// Let's look at the options: ip and port
	// A very basic usage of BOOST program options
	std::string varxmlfile,monxmlfile,sfavorites;
	// start adding options
	po::options_description desc("test_monitor options ");
	desc.add_options()
	        ("help,h", "An help message")
	        ("varxmlfile,f", po::value<std::string>(&varxmlfile)->default_value("HK-OSCVariables.xml"), "The variables definition XML filename")
	        ("monxmlfile,m", po::value<std::string>(&monxmlfile)->default_value("Monitor.xml"), "The monitor definition XML filename")
			("favname,s", po::value<std::string>(&sfavorites)->default_value("favorites"), "The name of the XML tree with the favorite definitions")
	;
	// Create the option map
	po::variables_map vm;
	// Parse command line for options
	po::store(po::parse_command_line(ac, av, desc), vm);
	// Store user options into the map and the local variables
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 0;
	}

	std::cout << " Using XML files: " << varxmlfile << ", " << monxmlfile << ", and favorites label: " << sfavorites << std::endl;
	 // Create empty property tree object
	pt::ptree tree;

	// Parse the XML into the property tree.
	pt::read_xml(varxmlfile, tree);
	std::cout << " Printing the XML file data \n";
	printtree(tree);
	boost::optional<pt::ptree &> my_tree=tree.get_child_optional(fcutils::monitor::xml::slowcontrol::VARIABLESROOT);
	boost::optional<pt::ptree &> valtree;
	if(my_tree){
		valtree=find_by_value_optional(*my_tree,"hkscvarHKEND");
	}
	std::cout << " Printing the hkscavarHKEND \n";
	if(valtree) printtree(*valtree);

	fcutils::monitor::value v;
	v.print();
	std::cout << "\n Fetching the data ...\n";
	fcutils::monitor::fetch_value(*valtree,v);
	v.print();
//	v.conversion_factor = 100.;
	std::cout << fcutils::terminal::ansi::color::BG_DEFAULT << " test terminal background ... " << fcutils::terminal::ansi::format::TXT_RESET << std::endl;
	double temp;
	while(std::cin>>temp){
		v=temp;
		v.print();
		// Just a test ... change the input data type
		auto intype=v.get_intype();
		v.set_intype(fcutils::monitor::inouttype::DOUBLE);
		v=temp;
		v.print();
		v.set_intype(intype);
		std::cout << v << std::endl;
	}
	fcutils::monitor::screen<fcutils::monitor::value, unsigned char> scrn;
//	pt::ptree monitor_tree;
//	pt::read_xml(monxmlfile,monitor_tree);
//	printtree(monitor_tree);
	scrn.fetch_all(monxmlfile);
	scrn.print_all();
	fcutils::monitor::adapt_terminal(scrn);
//	fcutils::monitor::list_values(scrn);
//	fcutils::monitor::list_labels(scrn);
	std::cout << scrn << std::endl;
	while(std::cin>>temp){
		fcutils::monitor::adapt_terminal(scrn);
		std::cout << scrn << std::endl;
	}
	return 0;
}

