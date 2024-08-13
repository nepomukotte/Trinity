///
///
/// @brief A new file
///
/// @file guarded_thread.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 20 ott 2021 
///

#ifndef GUARDED_THREAD_HPP_
#define GUARDED_THREAD_HPP_

/// Thanks to B.Stroustrup for the idea
#include <iostream>
#include <thread>

namespace fcutils{
namespace concurrency{

enum class THREAD_DTOR_POLICY {JOIN,DETACH};

template <bool Debug=false, THREAD_DTOR_POLICY DP=THREAD_DTOR_POLICY::JOIN> struct guarded_thread_base : std::thread {
	using std::thread::thread;
	using std::thread::operator=;

	guarded_thread_base() = default;
	// Copy ...
	guarded_thread_base(const std::thread & t) = delete;
	guarded_thread_base(const guarded_thread_base & t) = delete;
	// ... move
	guarded_thread_base(guarded_thread_base && t) noexcept =default;
	guarded_thread_base & operator=(guarded_thread_base && t) noexcept = default;
// 	Enable moving a standard thread into a guarded_thread
	guarded_thread_base(std::thread && t) noexcept: std::thread{std::move(t)}{};
	guarded_thread_base & operator=(std::thread && t) noexcept{
		std::thread::operator=(std::move(t));
		return *this;
	}

	constexpr bool is_debug() const { return Debug; }
	~guarded_thread_base() { if (joinable()) {
		// Choose the destruction policy
		if constexpr (DP==THREAD_DTOR_POLICY::JOIN){
			if constexpr (Debug) { std::cout<< "... joining ...\n"; }
			join();
		} else  {
			if constexpr (Debug) { std::cout<< "... detaching ...\n"; }
			detach();
		}
	}
	}
};

using guarded_thread=guarded_thread_base<false>;
using guarded_thread_debug=guarded_thread_base<true>;
// Joinable threads
using joinable_guarded_thread=guarded_thread_base<false>;
using joinable_guarded_thread_debug=guarded_thread_base<true>;
// Detachable threads
using detachable_guarded_thread=guarded_thread_base<false,THREAD_DTOR_POLICY::DETACH>;
using detachable_guarded_thread_debug=guarded_thread_base<true,THREAD_DTOR_POLICY::DETACH>;
// Instantiate the nodebug template here
template struct guarded_thread_base<false>;
template struct guarded_thread_base<false,THREAD_DTOR_POLICY::DETACH>;

} // namespace fcutils
} // namespace concurrency

#endif // GUARDED_THREAD_HPP_ 
