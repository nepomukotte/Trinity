///
///
/// @brief net_duplicator_manager, Map ostream duplicators according to their address, port and transport.
/// A key is: "ip_address:port transport". Transport is just a char: 'U' or 'T'.
/// For example  an ostream duplicators on the ip address 192.168.1.33, on port 514, transport TCP, will be represented by the key: "192.168.1.33:513 T".
/// The key case is always upper.
///
/// @file net_duplicator_manager.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 24 mar 2022 
///

#ifndef NET_DUPLICATOR_MANAGER_HPP_
#define NET_DUPLICATOR_MANAGER_HPP_
#include <memory>
#include <map>
#include <string>
#include <thread>

#include "fcutils/iostream/net_duplicator.hpp"
#include "fcutils/concurrency/lockable_container.hpp"

namespace fcutils{
namespace iostream{

template<typename BaseDuplicator=fcutils::iostream::base_net_duplicator,
		typename UdpDuplicator=fcutils::iostream::ostream_udp_duplicator,
		typename TcpDuplicator=fcutils::iostream::ostream_tcp_duplicator>
class net_duplicator_manager {
public:
	using base_net_duplicator_type	= BaseDuplicator;
	using tcp_net_duplicator_type	= TcpDuplicator;
	using udp_net_duplicator_type	= UdpDuplicator;
	using duplicator_ptr_type		= std::shared_ptr<base_net_duplicator_type>;
	using net_duplicator_map_type	= fcutils::concurrency::lockable_container<std::map,std::string,duplicator_ptr_type>;

	// typedef std::map<std::string,fcutils::iostream::base_net_duplicator *> udp_duplicator_map;
//	~net_duplicator_manager(){
//		std::for_each(udm_.begin(),udm_.end(), deleter<udp_duplicator_map::value_type>());
//	}

	bool enable(const char *c, std::ostream & o){
		return enable(std::string(c),o);
	}

	bool enable(std::string temp, std::ostream & o){
		bool done=false;
		// A string transformation is needed to be sure we are using the same case
		std::transform( temp.begin(), temp.end(), temp.begin(), tolower);
		std::stringstream ss(temp);
		// Lock the map to avoid data races
		std::lock_guard< net_duplicator_map_type > lck{duplicators_};
		auto it=duplicators_.find(ss.str());
		// If the duplicator has been already inserted, just enable it. Otherwise insert it
		// If is a TCP duplicator, a test of the connection is done before creating it.
		if(it==duplicators_.end()){
			std::string address;
			int port;
			std::getline(ss,address,':');
			ss>>port;
			char transport;
			ss >> transport;
	//		std::cout <<" New address: "<< address <<", and port: "<<port <<", transport "<< transport <<std::endl;
			duplicator_ptr_type duplicator=nullptr;
			switch(transport){
			case 'u':
				(duplicators_[ss.str()]=std::static_pointer_cast<base_net_duplicator_type>(std::make_shared<udp_net_duplicator_type>(address.c_str(),port)) )->enable(o);
				break;
			case 't':
				(duplicators_[ss.str()]=std::static_pointer_cast<base_net_duplicator_type>(std::make_shared<tcp_net_duplicator_type>(address.c_str(),port)) )->enable(o);
				break;
			default:
				break;
			}
		} else {
			it->second->enable(o);
			done=true;
		} // it==duplicators_.end()
		return done;
	}

	bool disable(const char *c, std::ostream & o){
		return disable(std::string{c},o);
	}

	bool disable( std::string temp, std::ostream &o){
		bool done=false;
		// A string transformation is needed to be sure we are using the same case
		std::transform( temp.begin(), temp.end(), temp.begin(), tolower);
		std::string ss(temp);
		auto it=duplicators_.find(ss);
//		std::cout << " Disabling UDP duplicator on: " << ss << std::endl;
		if(it!=duplicators_.end()){
//			std::cout << " Disabled " << std::endl;
			it->second->disable(o);
			done=true;
		}
		return done;
	}

	void disable_all(std::ostream &o){
		// Lock the map to avoid data races
		std::lock_guard< net_duplicator_map_type > lck{duplicators_};
		for( auto & [key,ptr]: duplicators_){
			ptr->disable(o);
		}
	}

	void enable_all(std::ostream &o){
		// Lock the map to avoid data races
		std::lock_guard< net_duplicator_map_type > lck{duplicators_};
		for( auto & [key,ptr]: duplicators_){
			ptr->enable(o);
		}
	}

private:
	net_duplicator_map_type duplicators_;

//	template<class Pair, class T> class apply_func{
//	public:
//		typedef void (T::*func)(std::ostream &);
//		apply_func(std::ostream & o, func f ): o_(o), f_(f){}
//		void operator()( Pair &p){
//			(p.second->*f_)(o_);
//		}
//	private:
//		func f_;
//		std::ostream &o_;
//	};
};

} // namespace iostream

} // namespace fcutils



#endif // NET_DUPLICATOR_MANAGER_HPP_ 
