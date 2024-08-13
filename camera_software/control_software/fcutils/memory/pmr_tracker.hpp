///
///
/// @brief To track allocation and deallocation of a PMR memory resource. Thanks to N. Josuttis for the idea and code.
///
/// @file pmr_tracker.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 7 dic 2021 
///

#ifndef PMR_TRACKER_HPP_
#define PMR_TRACKER_HPP_

#include <iostream>
#include <string>
#include <memory_resource>

namespace fcutils{
namespace memory{

class pmr_tracker : public std::pmr::memory_resource
{
 public:
  // we wrap the passed or default resource:
  explicit pmr_tracker(std::pmr::memory_resource* us
                    = std::pmr::get_default_resource())
   : upstream_{us} {
  }

  explicit pmr_tracker(std::string p,
                   std::pmr::memory_resource* us
                    = std::pmr::get_default_resource())
   : prefix_{std::move(p)}, upstream_{us} {
  }

 private:

   void* do_allocate(size_t bytes, size_t alignment) override {
    std::cout << prefix_ << "allocate " << bytes << " Bytes\n";
    void* ret = upstream_->allocate(bytes, alignment);
    return ret;
  }

  void do_deallocate(void* ptr, size_t bytes, size_t alignment) override {
    std::cout << prefix_ << "deallocate " << bytes << " Bytes\n";
    upstream_->deallocate(ptr, bytes, alignment);
  }

  bool do_is_equal(const std::pmr::memory_resource& other) const noexcept
   override {
    // same object?:
    if (this == &other) return true;
    // same type and prefix and equal upstream?:
    auto op = dynamic_cast<const pmr_tracker*>(&other);
    return op != nullptr && op->prefix_ == prefix_
            && upstream_->is_equal(other);
  }

  std::pmr::memory_resource* upstream_;  // wrapped memory resource
  std::string prefix_{};

};


} // namespace memory
} // namespace fcutils




#endif // PMR_TRACKER_HPP_ 
