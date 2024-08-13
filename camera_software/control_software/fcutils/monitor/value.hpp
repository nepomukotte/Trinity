///
///
/// @brief A new file
///
/// @file value.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 26 mar 2022 
///

#ifndef VALUE_HPP_
#define VALUE_HPP_
#include <iostream>
#include <string>
#include <variant>
#include <optional>
#include <type_traits>
#include <utility>

#include <boost/property_tree/ptree.hpp>
#include <boost/format.hpp>

#include "fcutils/terminal/ansi.hpp"
#include "fcutils/type_traits/type_traits.hpp"
#include "fcutils/string/change_case.hpp"
#include "fcutils/monitor/string_to_object.hpp"
#include "fcutils/monitor/color_map.hpp"
#include "fcutils/monitor/endianess_map.hpp"
#include "fcutils/monitor/inouttype_map.hpp"
#include "fcutils/monitor/outformat_map.hpp"

namespace fcutils{
namespace monitor{

template< class Value> Value mask_and_shift(Value v, std::optional<std::size_t> mask=std::nullopt,
        std::optional<std::size_t> rshift=std::nullopt){
        if constexpr (std::is_integral<Value>::value){
                if(mask) v= v& *mask;
                if(rshift) v= v >> *rshift;
        }
        return v;
};


struct converter{
	converter(double f=1.,double off=0., std::size_t bm=~0, std::size_t bs=0): factor(f), offset(off), bytemask(bm), bitshift(bs) {}
	double factor;
	double offset;
	std::size_t bytemask;
	std::size_t bitshift;
	template <class Return, class Value > Return operator()( Value v){
		v=mask_and_shift(v,bytemask,bitshift);
		return static_cast<Return>(v)*static_cast<Return>(factor)+static_cast<Return>(offset);
	}
	template <class Return, class Value > Return convert( Value v){
		return operator()(v);
	}
};



struct value{
	using tree_type 			= boost::property_tree::ptree;
	using self_type 			= value;
	using value_type			= fcutils::monitor::value_type; // from inouttype_map.hpp
	using converted_value_type	= fcutils::monitor::value_type; // from inouttype_map.hpp
	using color_type			= fcutils::terminal::ansi::color;
	using typelist_type 		= fcutils::monitor::typelist_type; // from inouttype_map.hpp
	using outtype_type			= fcutils::monitor::inouttype;
	using intype_type			= fcutils::monitor::inouttype;
	using endianess_type		= fcutils::monitor::endianess;
	using outformat_type		= std::ios_base::fmtflags;

	std::string name;
	std::string id;
	std::string label;
	size_t offset{0};
	size_t size{0};
	endianess_type endian{endianess_type::BIGENDIAN};
	outformat_type outfmt{std::ios_base::dec};
	std::optional<color_type> alarmcol;
	std::optional<converted_value_type> max;
	std::optional<converted_value_type> min;
	std::optional<double> conversion_factor;
	std::optional<double> conversion_offset;
	std::optional<std::size_t> conversion_bytemask;
	std::optional<std::size_t> conversion_bitshift;
	bool alarm{false};
private:
	value_type val_{0};
	converted_value_type converted_val_{0.};
	outtype_type outtype_{outtype_type::DOUBLE};
	intype_type intype_{intype_type::INT};

public:

	value() = default;
	value( const value &) = default;
	value( value &&) = default;
	value & operator=( const value & ) = default;
	value & operator=( value && ) = default;
	template<typename Val> self_type & operator=(const Val v){
		// Search for the type in the list. If present then assign and convert
		if constexpr (fcutils::type_traits::list_search<Val,typelist_type>::index >0 ){
			// The type is in the list. Assign it
			std::visit([&v]( auto & val_ ) mutable { val_=static_cast<fcutils::type_traits::remove_cvref_t<decltype(val_)>>(v);}, val_);
			// Enable the alternative according to the output type
			enable_alternative(outtype_,converted_val_);
			// Convert it if both mask and bitshift are set. Check is done in mask_and_shift
			auto converted=mask_and_shift(v,conversion_bytemask,conversion_bitshift);
			// Convert if conversion factor and offset are set
			std::visit(
					[converted,conv_f=conversion_factor,conv_o=conversion_offset](auto & v){
						using value_type = fcutils::type_traits::remove_cvref_t<decltype(v)>;
						v = (conv_f)?static_cast<value_type>(converted)*static_cast<value_type>(*conv_f):static_cast<value_type>(converted);
						v = (conv_o)?v+static_cast<value_type>(*conv_o):v;
			}, converted_val_);
			// Check for limits if set
			alarm=(max)?(converted_val_>*max):false;
			alarm=(min)?(alarm||(converted_val_<*min)):alarm;
		}
		return *this;
	}

	auto get_value() const {
		return val_;
	}

	auto get_converted_value() const {
		return converted_val_;
	}

	void set_intype(intype_type it){
		intype_=it;
		// Reset the alternative hold in the variant
		enable_alternative(it,val_);
	}

	void set_outtype(outtype_type ot){
		outtype_=ot;
		// Reset the alternative hold in the variant
		enable_alternative(ot,converted_val_);
		if(max) {
			double temp;
			std::visit([&temp](auto & v)mutable{ temp=static_cast<fcutils::type_traits::remove_cvref_t<decltype(temp)>>(v);},*max);
			set_alternative(ot,*max,temp);
		}
		if(min) {
			double temp;
			std::visit([&temp](auto & v)mutable{ temp=static_cast<fcutils::type_traits::remove_cvref_t<decltype(temp)>>(v);},*min);
			set_alternative(ot,*min,temp);
		}
	}

	void set_inouttype(intype_type it, outtype_type ot){
		set_intype(it);
		set_outtype(ot);
	}

	template <typename Data> void set_value( const Data d){
		operator= (d);
	}

	template <typename Data> void set_value( const Data *d){
//		operator= (d);
		std::cout << " In the set_value with pointer ... to be implemented! "<< std::endl;
	}

	//
	auto get_intype() const { return intype_;}
	auto get_outtype() const { return outtype_;}

	std::ostream & to_stream(std::ostream &os) const {
		auto m=fcutils::terminal::ansi::color::BG_DEFAULT;
		if(alarmcol) m=(alarm)?*alarmcol:m;
		os << m;
		std::visit([&os]( auto & v ) { os << v; },converted_val_);
		os << fcutils::terminal::ansi::format::TXT_RESET;
		return os;
	}

    boost::format & to_format(boost::format & fmt) const {
		std::visit([&fmt]( auto & v ) mutable { fmt % v; },converted_val_);
    	return fmt;
    }

    boost::format & to_format(boost::format & fmt, std::size_t width, std::size_t precision) const {
    		std::stringstream ss;
    		// First set the background color if an alarm occurred
    		auto m=fcutils::terminal::ansi::color::BG_DEFAULT;
            if(alarmcol) m=(alarm)?*alarmcol:m;
            ss << m;
            // Then format according to the output type
    		switch(outtype_){
    		case outtype_type::DOUBLE:
    		case outtype_type::FLOAT:
    			// floating point case, build the format string with width and precision: %w.pf
    			ss<< '%' << width << '.' << precision << 'f';
    			break;
    		case outtype_type::INT:
    			// Integer type, check if the hex has been request
    		{
    			char f=(outfmt==std::ios_base::hex)?'X':'i';
    			ss<< "%#" << width << f;
    			break;
    		}
    		default:
    			// all the others, set just the width
    			ss<< '%' << width << 'd';
    			break;
    		}
    		// then set back the background color
    		ss << fcutils::terminal::ansi::color::BG_DEFAULT;
    		// Now copy the format string into the provided format object
    		fmt=boost::format(ss.str().c_str());
    		// Save the value into the format object and return it
    		std::visit([&fmt]( auto & v ) mutable { fmt % v; },converted_val_);
    		return fmt;
    }

    void print(std::ostream & os=std::cout) const {
		os << " Value name: " << name <<", id: "<< id <<", label: " << label <<'\n';
		os << " Memory offset: " << offset <<", size (in Byte) :" << size << ", intype: " <<
				fcutils::type_traits::enum_to_underlying(intype_) << ", outtype: " << fcutils::type_traits::enum_to_underlying(outtype_) << '\n';
		if(conversion_bytemask) {
			os << " Bytemask: " << std::hex << std::showbase << *conversion_bytemask << std::noshowbase << std::dec;
		} else {
			os << " No Bytemask ";
		}
		if(conversion_bitshift){
			os << ", bitshift: "<< *conversion_bitshift << '\n';
		} else {
			os << ", no bitshift. \n";
		}
		if(conversion_factor){
			os << " Conversion factor: "<< *conversion_factor;
		} else {
			os << " No conversion factor, ";
		}
		if(conversion_offset){
			os << ", conversion offset: " << *conversion_offset << '\n';
		} else {
			os << " No conversion offset. \n";
		}
		os << " Endianness: " << fcutils::type_traits::enum_to_underlying(endian) << '\n'; //<< ", out format: " << outfmt << '\n';
		std::string temp{" Raw data value: "};
		auto visitor=[&temp,&os](auto &v){ os << temp << v; };
		std::visit(visitor,val_);
		temp= ". Converted data value: ";
		std::visit(visitor,converted_val_);
		if(min){
			temp= "\n Min: ";
			std::visit(visitor,*min);
		} else {
			os << "No Min value set. ";
		}
		if(max){
			temp= ". Max: ";
			std::visit(visitor,*max);
		} else {
			os << " No Max value set. \n";
		}
		std::cout << std::boolalpha << ". Alarm: " << alarm << std::noboolalpha << std::endl;
	}

	template <typename Tree> bool fetch(Tree && pt){
		return fetch_value(std::forward<Tree>(pt),*this);
	}
};

std::ostream & operator <<(std::ostream &os, const value & v){
	// Backup the actual flags
	std::ios_base::fmtflags ff{os.flags()};
	// Set the format based on the value one
	os.setf (v.outfmt , std::ios::basefield);
	// Print the value
	v.to_stream(os);
	os.flags(ff);
	return os;
}


bool fetch_value( value::tree_type & pt, value & v){
		const std::size_t FIELDWIDTH{6};
		bool success=false;
		v.name=pt.get_child(xml::slowcontrol::variable::field::NAMEOPT).get_value("noname");
		v.id=pt.get_child(xml::slowcontrol::variable::field::IDOPT).get_value("noid");
		v.label=pt.get_child(xml::slowcontrol::variable::field::LABELOPT).get_value("nolabel");
		if(v.label.size()<FIELDWIDTH) v.label=std::string(FIELDWIDTH-v.label.size(),' ') + v.label;
		v.offset=pt.get_child(xml::slowcontrol::variable::field::OFFSETOPT).get_value(-1);
		v.size=pt.get_child(xml::slowcontrol::variable::field::SIZEOPT).get_value(0);
		std::string temp=pt.get_child(xml::slowcontrol::variable::field::INTYPEOPT).get_value(xml::slowcontrol::variable::value::SDOUBLE);
		// default in and out type to int
		v.set_inouttype(inouttype::INT, inouttype::INT);
		auto ittemp=string_to_object(temp,STRING_TO_INOUTTYPE);
		v.set_intype(ittemp);
		auto ptvar=pt.get_child_optional(xml::slowcontrol::variable::field::OUTTYPEOPT);
		if(ptvar) {
			temp=ptvar->get_value(xml::slowcontrol::variable::value::SDOUBLE);
			auto ottemp=string_to_object(temp,STRING_TO_INOUTTYPE);
			v.set_outtype(ottemp);
		}
		ptvar.reset();
		temp=pt.get_child(xml::slowcontrol::variable::field::ENDIANNESSOPT).get_value(xml::slowcontrol::variable::value::SBIGENDIAN);
		v.endian= string_to_object(temp, STRING_TO_ENDIANESS);
		temp=pt.get_child(xml::slowcontrol::variable::field::FORMATOPT).get_value(xml::slowcontrol::variable::value::SDEC);
		v.outfmt= string_to_object(temp, STRING_TO_OUTFORMAT);
		ptvar=pt.get_child_optional(xml::slowcontrol::variable::field::ALARMCOLOROPT);
		if(ptvar) {
			temp=ptvar->get_value(xml::slowcontrol::variable::value::SRED);
			v.alarmcol=string_to_object(temp, STRING_TO_BG_COLOR);
		}
		ptvar.reset();
		ptvar=pt.get_child_optional(xml::slowcontrol::variable::field::MAXOPT);
		if(ptvar) {
			auto max=ptvar->get_value<double>();
			// Enable the optional
			v.max=max;
			set_alternative(v.get_outtype(),*v.max,max);
		}
		ptvar.reset();
		ptvar=pt.get_child_optional(xml::slowcontrol::variable::field::MINOPT);
		if(ptvar) {
			auto min=ptvar->get_value<double>();
			// Enable the optional
			v.min=min;
			set_alternative(v.get_outtype(),*v.min,min);
		}
		ptvar.reset();
		ptvar=pt.get_child_optional(xml::slowcontrol::variable::field::CONVERSIONFACTOROPT);
		if(ptvar) {
			v.conversion_factor=ptvar->get_value<double>();
		}
		ptvar.reset();
		ptvar=pt.get_child_optional(xml::slowcontrol::variable::field::CONVERSIONOFFSETOPT);
		if(ptvar) {
			v.conversion_offset=ptvar->get_value<double>();
		}
		ptvar.reset();
		ptvar=pt.get_child_optional(xml::slowcontrol::variable::field::CONVERSIONBYTEMASKOPT);
		if(ptvar) {
			std::string temp=ptvar->get_value(xml::slowcontrol::variable::value::DEFBYTEMASK);
			std::stringstream ss;
			ss << std::hex << temp;
			ss >> *v.conversion_bytemask;
		}
		ptvar.reset();
		ptvar=pt.get_child_optional(xml::slowcontrol::variable::field::CONVERSIONBITSHIFTOPT);
		if(ptvar){
			v.conversion_bitshift=ptvar->get_value<std::size_t>();
		}

		success=true;
		return success;
	}

} // namespace monitor
} // namespace fcutils



#endif // VALUE_HPP_ 
