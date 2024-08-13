///
///
/// @brief An helper function to print a variable numbers of values using a delimiter, an header and trailer
///
/// @file tuple_to_string.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 3 feb 2022 
///

#ifndef UTILITY_PRINT_VALS_HPP_
#define UTILITY_PRINT_VALS_HPP_
#include <iosfwd>

namespace fcutils{

template<char Delimiter=',', char Header=' ', char Trailer='\n', typename ... Arg>
inline void print_vals(std::ostream & o, Arg ... args){
		auto n= sizeof ...(args);
		auto print=[&o](const auto & v, auto count, const auto & delimiter){ o<<v; if(count>0) o<<delimiter<<' '; };
        o << Header;
		(print(args,--n,Delimiter), ...);
        o << Trailer;
}

} // namespace fcutils

#endif // UTILITY_PRINT_VALS_HPP_ 
