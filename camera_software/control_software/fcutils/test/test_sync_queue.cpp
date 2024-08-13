///
///
/// @brief A new file
///
/// @file test_sync_queue.cpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 8 nov 2021 
///

#include <concurrency/guarded_thread.hpp>
#include <concurrency/sync_queue.hpp>
#include <string>
#include <iostream>
#include <utility>


using my_queue=fcutils::concurrency::sync_queue<std::string,false>;

struct consumer{
	my_queue & head, & user;
	int i{0};
	consumer(my_queue & q, my_queue & u) :head{q}, user{u} {};
	void operator()(){
		++i;
		std::string msg;
		head.get(msg);
		if(msg[0]==eotx) {
			std::cout << " Consumer end of transmission detected ... \n";
			user.put(std::move(msg));
			return;
		}
		std::cout << " Consumer i: "<< i<<", msg: "<<msg <<std::endl;
		(*this)();
	};
	const char eotx{0x04};
};

struct producer{
	my_queue & tail;
	producer(my_queue & q) :tail(q) {};
	void operator()(){
		std::string msg;
		std::cout << " Producer. Enter a message: \n";
		while(std::getline(std::cin,msg)){
			std::cout << " Producer. Pushing message: "<< msg << '\n';
			tail.put(std::move(msg));
			std::cout << " Producer. Enter a message: \n";
		}
		msg=eotx;
		tail.put(std::move(msg));
		return;
	};
	const char eotx{0x04};
};

int main(){
	my_queue mq,uq;

	consumer c{mq,uq};
	producer p{mq};
	fcutils::concurrency::guarded_thread tc{c};
	fcutils::concurrency::guarded_thread tp{p};
	std::string msg;
	uq.get(msg);
	std::cout << " A stop message was sent ... stopping the main... " << std::endl;
	return 0;
}

