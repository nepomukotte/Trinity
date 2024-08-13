///
///
/// @brief The finally pattern. An alternative to exception throw and RAII
///
/// @file finally.hpp
/// @author Francesco S. Cafagna (francesco.cafagna@ba.infn.it)
/// @date 16 dic 2021 
///

#ifndef UTILITY_FINALLY_HPP_
#define UTILITY_FINALLY_HPP_

namespace fcutils{

template <typename Funct> struct final_action{
        final_action(Funct f): clean{f} {};
        ~final_action(){ clean(); };
        Funct clean;
};

template<typename Funct> final_action<Funct> finally(Funct f){
        return final_action<Funct>(f);
}

} // namespace fcutils

#endif // UTILITY_FINALLY_HPP_ 
