#ifndef UTILS_DATE_STDTMHELPERS_H_INCLUDED
#define UTILS_DATE_STDTMHELPERS_H_INCLUDED

#include <istream>
#include <ostream>
#include <ctime>

bool operator==(const std::tm& t1, const std::tm& t2);
std::istream& operator>>(std::istream& is, std::tm& dt);
std::ostream& operator<<(std::ostream& os, const std::tm& dt);

#endif // UTILS_DATE_STDTMHELPERS_H_INCLUDED
