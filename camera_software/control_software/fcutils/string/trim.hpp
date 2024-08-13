///
///
/// @brief A quick and dirty trim function using C++17 string_view
///
/// @file trim.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 15 mar 2022 
///

#ifndef TRIM_HPP_
#define TRIM_HPP_

#include <string_view>
#include <algorithm>

namespace fcutils{
namespace string{

std::string_view trim(std::string_view s)
{
    s.remove_prefix(std::min(s.find_first_not_of(" \t\r\v\n"), s.size()));
    s.remove_suffix(std::min(s.size() - s.find_last_not_of(" \t\r\v\n") - 1, s.size()));

    return s;
}

} // namespace string
} // namespace fcutils




#endif // TRIM_HPP_ 
