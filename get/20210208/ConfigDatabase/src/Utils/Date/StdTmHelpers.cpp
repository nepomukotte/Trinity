#include "StdTmHelpers.h"
#include "../Stream/Require.h"
#include <iomanip>
#include <cstring>

/**
 * Compares 2 timestamps for equality.
 */
bool operator==(const std::tm& t1, const std::tm& t2)
{
	return t1.tm_sec == t2.tm_sec and
		t1.tm_min == t2.tm_min and
		t1.tm_hour == t2.tm_hour and
		t1.tm_mday == t2.tm_mday and
		t1.tm_mon == t2.tm_mon and
		t1.tm_year == t2.tm_year;
}

/**
 * Extracts a timestamp using the ISO 8601 format from a stream.
 * @note only UTC times are supported
 */
std::istream& operator>>(std::istream& is, std::tm& dt)
{
	memset(&dt, 0, sizeof dt);
	is	>> dt.tm_year >> require('-')
		>> dt.tm_mon >> require('-')
		>> dt.tm_mday
		>> require('T')
		>> dt.tm_hour >> require(':')
		>> dt.tm_min >> require(':')
		>> dt.tm_sec
		>> require('Z');
	dt.tm_year -= 1900;
	dt.tm_mon -= 1;
	return is;
}

/**
 * Inserts a timestamp using the ISO 8601 format on a stream.
 * @note only UTC times are supported
 */
std::ostream& operator<<(std::ostream& os, const std::tm& dt)
{
	char oldFill = os.fill('0');
	os  << std::setw(4) << dt.tm_year + 1900 << '-'
		<< std::setw(2) << dt.tm_mon + 1 << '-'
		<< std::setw(2) << dt.tm_mday
		<< 'T'
		<< std::setw(2) << dt.tm_hour << ':'
		<< std::setw(2) << dt.tm_min << ':'
		<< std::setw(2) << dt.tm_sec
		<< 'Z';
	os.fill(oldFill); // restores the fill
	return os;
}
