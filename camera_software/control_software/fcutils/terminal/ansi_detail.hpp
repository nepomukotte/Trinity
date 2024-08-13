///
///
/// @brief A new file
///
/// @file ansi_details.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 1 feb 2022 
///

#ifndef ANSI_DETAIL_HPP_
#define ANSI_DETAIL_HPP_
#include <iosfwd>
#include <tuple>

namespace fcutils{
namespace terminal{
namespace ansi{

using trailer_type=char;
using parameter_type=int;

enum class csi: trailer_type {
	ESCAPE='\033', // Needed for the CSI sequence "ESC[ # CURSOR_COMMAND SGR", being # an optional value
	SQUAREBRACKET='[',
	COLON=':',
	SGR='m' // Select Graphic Rendition
};

// Stream of the csi class
std::ostream & operator<<(std::ostream & o, const csi & c){
	return o << static_cast<trailer_type>(c);
}

/// A simple class to map an ANSI parameter to a type
template<typename C, C Parameter, bool SGR=true> struct par2type{
	enum : C { PARAMETER=Parameter };
};
/// The case with the SGR at the end of the CSI
template<typename C=trailer_type, C Parameter, template <typename, C,bool> class T>
std::ostream & operator<<(std::ostream & o, const T<C,Parameter,true> & cmd){
	return o << csi::ESCAPE <<
			csi::SQUAREBRACKET <<
			static_cast<C>(cmd.PARAMETER) <<
			csi::SGR;
}
// No SGR at the end of the CSI
template<typename C=trailer_type, C Parameter, template <typename, C,bool> class T>
std::ostream & operator<<(std::ostream & o, const T<C,Parameter,false> & cmd){
	return o << csi::ESCAPE <<
			csi::SQUAREBRACKET <<
			static_cast<C>(cmd.PARAMETER) ;
}

template<typename Code, Code c > struct escape_sequence{
	template <typename ... Arg> decltype(auto) operator()(std::ostream & o, Arg ... args){
		auto n= sizeof ...(args);
		auto print=[&o](const auto & v, auto count, const auto & delimiter){ o<<v; if(count>0) o<<delimiter; };
        o << csi::ESCAPE << csi::SQUAREBRACKET;
		(print(args,--n,';'), ...);
        o << c;
	}
};


} // namesapce ansi

} // namesapace terminal

} // namespace fcutils



#endif // ANSI_DETAIL_HPP_ 
