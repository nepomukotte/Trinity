///
///
/// @brief A new file
///
/// @file return_policies.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 1 nov 2021 
///

#ifndef RETURN_POLICIES_HPP_
#define RETURN_POLICIES_HPP_

#include <utility> // for std::forward

namespace fcutils {
namespace checker {

/// Simply return the value of the checker as it is
template<typename Checker, typename Return> class do_nothing
{
public:
  using return_type=Return;
  using checker_type=Checker;
  template<typename Data> return_type operator()( const checker_type & r, const Data *, std::size_t)
  {
    return r;
  }
};

/// In case of bool and size_t, the number of word to be read next is returned
/// Just a specialization of the general case
template<> class do_nothing<bool,std::size_t>
{
public:
  using return_type=std::size_t;
  using checker_type=bool;
  template<typename Data> return_type operator()( const checker_type & r, const Data *, std::size_t)
  {
	  // If the Checker condition is true, return zero word to be read, otherwise return 1
    return r?0:1;
  }
};

/// Return the number of words to be read, decoding a word of an header
template<typename Checker, typename Return> class decode_header;
template<> class decode_header<bool,std::size_t>
{
public:
  using self_type=decode_header<bool,std::size_t>;
  using return_type=std::size_t;
  using checker_type=bool;
  decode_header(std::size_t h, std::size_t o):header_size_{h}, offset_{o}, decoded_{} {};
  decode_header(const self_type & d) = default;
  self_type & operator=(const self_type & d) = default;

  template<typename Data> return_type operator()( const checker_type & r, Data *d, std::size_t l)
  {
	  // If the number of words read is the offset then return the data casted to the return_type
	  if(l==offset_) return decoded_=static_cast<return_type>(*d );
	  // If the number of words is less equal than the header_size, but is not the offset, then return 1
	  if(l<=header_size_) return 1;
	  // If the number of words is larger than the header size then stop reading
	  return 0;
  }
private:
  std::size_t header_size_,offset_;
  return_type decoded_;
};

/// A function to return a policy object, useful to pass a policy object that needs initialization
template<typename Policy, typename ... Arg> inline Policy return_policy(Arg && ... args){
	return Policy{std::forward(args)...};
}

} // namespace checker
} // namespace fcutils




#endif // RETURN_POLICIES_HPP_ 
