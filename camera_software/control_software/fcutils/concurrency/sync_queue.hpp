///
///
/// @brief A new file
/// Thanks to Bjarne Stroustrup. The C++ Programming Language (4th Edition)
///
/// @file sync_queue.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 8 nov 2021 
///

#ifndef SYNC_QUEUE_HPP_
#define SYNC_QUEUE_HPP_

#include <mutex>
#include <condition_variable>
#include <chrono>
#include <list>
//#include <system_error>
#include <utility>

namespace fcutils{
namespace concurrency{

template<typename T, bool NotifyAll=false > class sync_queue {
public:
	using message_type=T;
	using list_type=std::list<message_type>;

	void put(const message_type& val){
		std::lock_guard<std::mutex> lck(mtx_);
		q_.push_back(val);
		if constexpr (NotifyAll) {
			cond_.notify_all();
		} else {
			cond_.notify_one();
		}
	}

	void put(message_type&& val){
		std::lock_guard<std::mutex> lck(mtx_);
		q_.push_back(std::move(val));
		if constexpr (NotifyAll) {
			cond_.notify_all();
		} else {
			cond_.notify_one();
		}
	}

	void get(message_type& val){
		std::unique_lock<std::mutex> lck(mtx_);
		cond_.wait(lck,[this]{ return !q_.empty(); });
		val=q_.front();
		q_.pop_front();
	}

//	void get(T& val, std::chrono::steady_clock::duration d) {
//		std::unique_lock<std::mutex> lck(mtx_);
//		bool not_empty = cond_.wait_for(lck,d,[this]{ return !q_.empty(); });
//		if (not_empty) {
//			val=q_.front();
//			q_.pop_front();
//		} else throw std::system_error{std::generic_category().default_error_condition(ETIMEDOUT),"Sync_queue: get() timeout"};
//	}

	bool get(message_type& val, std::chrono::steady_clock::duration d) {
		std::unique_lock<std::mutex> lck(mtx_);
		bool not_empty = cond_.wait_for(lck,d,[this]{ return !q_.empty(); });
		if (not_empty) {
			val=q_.front();
			q_.pop_front();
		}
		return not_empty;
	}

private:
	std::mutex mtx_;
	std::condition_variable cond_;
	list_type q_;

};

} // namespace concurrency
} // namespace fcutils

#endif // SYNC_QUEUE_HPP_ 
