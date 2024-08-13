/*
 * $Id$
 */

#include "get/rc/CoBoConfigParser.h"
#include <CCfg/Object.h>
#include <CCfg/CConfig.h>
#include <boost/test/auto_unit_test.hpp>
#include <iostream>

namespace get
{
namespace rc
{

BOOST_AUTO_TEST_CASE(AgetPeakingTimeIndex_test)
{
	// AGET existing peaking time
	size_t peakingTime = 117;
	uint16_t expectedIndex = 1;
	uint16_t actualIndex = CoBoConfigParser::peakingTimeIndex(peakingTime, 0x202);
	BOOST_CHECK_EQUAL(expectedIndex, actualIndex);

	// ASTRE undefined peaking time
	peakingTime = 117;
	expectedIndex = 6;
	actualIndex = CoBoConfigParser::peakingTimeIndex(peakingTime, 0x203);
	BOOST_CHECK_EQUAL(expectedIndex, actualIndex);

	// AGET undefined peaking time
	peakingTime = 3985;
	expectedIndex = 8;
	actualIndex = CoBoConfigParser::peakingTimeIndex(peakingTime, 0x202);
	BOOST_CHECK_EQUAL(expectedIndex, actualIndex);

	// ASTRE existing peaking time
	peakingTime = 3985;
	expectedIndex = 13;
	actualIndex = CoBoConfigParser::peakingTimeIndex(peakingTime, 0x203);
	BOOST_CHECK_EQUAL(expectedIndex, actualIndex);
}


} /* namespace rc */
} /* namespace get */

