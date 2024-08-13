///
///
/// @brief A new file
///
/// @file test_consumer.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 23 feb 2022 
///
#include <iostream>
#include <variant>
#include <map>
#include <type_traits>
#include <memory> // for the smart pointer facilities
#include <chrono>
#include <future>

#include "fcutils/networking/session.hpp" // for the pint and pint_byte consumer
#include "fcutils/utility/utility.hpp" // for the overload obkect
#include "fcutils/type_traits/type_traits.hpp"
#include "fcutils/concurrency/sync_queue.hpp"

namespace fcn=fcutils::networking;

using print_t=fcn::print<char>;
using print_byte_t=fcn::print_byte<char>;

template<typename T> struct is_shared_ptr : std::false_type {};
template<typename T> struct is_shared_ptr<std::shared_ptr<T>> : std::true_type {};

struct no_copy{
	no_copy(){};
	no_copy(const no_copy &) =delete;
	no_copy(no_copy &&) noexcept =default;
	no_copy & operator=(no_copy &&) noexcept = default;
	void operator()(const char *, std::size_t ){ std::cout << " Here I am: "<< std::hex << (void*)this << std::dec << std::endl; }
};

template < typename ... T> struct consumer_trait{
public:
	using self_type=consumer_trait<T ...>;
	using data_type=char;
//	template<typename ... Arg> consumer_trait(Arg && ... args): var{std::forward<Arg>(args) ... } {}

	consumer_trait(){};
	template <typename A> consumer_trait(A && a):var{std::forward<A>(a)}{};
//	consumer_trait(self_type && a):var{std::move(a.var)}{};
	self_type & operator=(self_type && a) noexcept = default;


	consumer_trait(self_type &&) = default;

	void consume(const data_type *d, std::size_t s){
		auto ovl=fcutils::overload{
//			[d,s](auto * a) { a->operator()(d,s); },
//			[d,s](auto & a) { a.operator()(d,s); },
			[d,s](auto & a) {
				if constexpr (is_shared_ptr<fcutils::type_traits::remove_cvref_t<decltype(a)>>::value) {
					a->operator()(d,s);
				} else {
					a(d,s); }
				}
		};
		std::visit(ovl,var);
	}

	bool operator()(const std::string & s){ consume(s.c_str(),s.size()); return true;}

	std::variant<T ... > var{};
};


using consume_call_type=bool (const std::string &);
using server_type=fcutils::concurrency::sync_queue<std::packaged_task<consume_call_type>>;

using namespace std::chrono_literals;

class tt{
public:
        tt(const std::chrono::seconds &s=5s): wait_{s}, this_id_{std::this_thread::get_id()}{};
        bool operator()(const std::string & s){
                std::this_thread::sleep_for(wait_);
                std::cout << " This_thread: " << this_id_
                << ", actual thread:" << std::this_thread::get_id()<<", string: " << s << std::endl;
                return true;
        }

//        ~t(){
//                std::future_status status;
//                do {
//                        switch(status = f_.wait_for(1s); status) {
//                                case std::future_status::deferred: std::cout << "deferred\n"; break;
//                                case std::future_status::timeout: std::cout << "timeout\n"; break;
//                                case std::future_status::ready: std::cout << "ready!\n"; break;
//                        }
//                 } while (status != std::future_status::ready);
//
//        }
//        void execute( const std::string & s){
//                f_=std::async(std::launch::async, &t::operator(),this,s);
//        }
private:
//        std::future<void> f_;
        std::thread::id this_id_;
        std::chrono::seconds wait_;
};


int main(){
	using variant_t=std::variant< std::shared_ptr<print_t>, std::shared_ptr<print_byte_t> >;
	variant_t var;
    std::string t{"Ciccillo2"};
    auto ovl=fcutils::overload{
    			[&t](auto a) { a->operator()(t.c_str(),t.size()); }
     };
    std::map<int,variant_t> consumers;
    var=std::make_shared<print_t>();
    consumers[0]=var;
//    std::get<0>(var)->operator()(t.c_str(),t.size());
//    std::visit(ovl,var);
    var=std::make_shared<print_byte_t>();
    consumers[1]=var;
    for(auto & [key,consumer]: consumers) std::visit(ovl,consumer);
    t="Ciccillo3";
    for(auto & [key,consumer]: consumers) std::visit(ovl,consumer);
//    std::visit(ovl,var);
//    std::get<1>(var)->operator()(t.c_str(),t.size());
    consumer_trait< std::shared_ptr<print_t>, std::shared_ptr<print_byte_t>, print_t, no_copy > ct;
    std::map<int, decltype(ct)::self_type > consumers2;
//    ct.var=std::make_shared<print_t>();
//    consumers2[0]=ct;
//    ct.var=std::make_shared<print_byte_t>();
//    consumers2[1]=ct;
//    ct.var=print_t{};
//    consumers2[2]=ct;
    no_copy junk;
    junk(t.c_str(),t.size());
//    ct.var=std::move(junk);
    consumers2.emplace(0,std::make_shared<print_t>());
    consumers2.emplace(1,std::make_shared<print_byte_t>());
    consumers2.emplace(2,print_t{});
    consumers2.emplace(3,std::move(junk));
    consumers2.emplace(4,std::make_shared<print_t>());
    t="Ciccillo4";
    for(auto & [key,consumer]: consumers2) consumer.consume(t.c_str(),t.size());

    server_type server;
//    for(auto & [key,consumer]: consumers2) {
    while(std::getline(std::cin,t)){
//    	server_type::message_type a(std::move(consumer));
    	std::packaged_task<bool(const std::string &)> a([](const std::string & s){std::cout << " Sono nel packaged_task ...: " << s << std::endl; return true;});
    	auto f=a.get_future();
    	a(t);
//    	tt junk;
//    	junk(t);
//    	server.put(std::move(a));
    	std::cout << std::boolalpha << f.get() << std::endl;
    }

	return 0;
}

