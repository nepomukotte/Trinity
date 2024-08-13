///
///
/// @brief A lockable container class
///
/// @file locable_container.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 10 mar 2022 
///

#ifndef LOCKABLE_CONTAINER_HPP_
#define LOCKABLE_CONTAINER_HPP_
#include <mutex>

namespace fcutils{
namespace concurrency{

// A lockable container class.
// It is just a container with a mutex added.
// Operation on mutex are forced through the public class interface.
// Just the container contructors are inherited
template <template <typename ... Arg> class Container, typename ... T>  class lockable_container: public Container<T ... >{
private:
        mutable std::mutex mtx_;
public:
        using container_type=Container<T ...>;
        // Inheriting constructors since C++11
        using container_type::Container;
        using container_type::operator=;
        void lock(){mtx_.lock();}
        void unlock(){ mtx_.unlock();}
        auto try_lock(){ return mtx_.try_lock();}
};

} // namespace concurrency
} // namespace fcutils



#endif // LOCKABLE_CONTAINER_HPP_ 
