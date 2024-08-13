/*
 * $Id: AgetLsbThresholdTest.cpp 1700 2017-09-15 12:53:54Z psizun $
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

BOOST_AUTO_TEST_CASE(Reg8Lsb_test)
{
	const std::string regName = "reg8";
	const uint64_t threshold = 0xF;
	// Check bitmask for individual channels
	for (size_t channelId=0; channelId < 68; channelId++)
	{
		// Create configuration
		CCfg::Object aget("Aget");
		CCfg::CConfig agetConfig(&aget);
		agetConfig("channel", channelId)("LSBThresholdValue").setIntValue(threshold);
		// Compute actual bitmask
		CCfg::CConfig::Iterator itChan = agetConfig.iterate("channel");
		itChan.next();
		uint64_t actualMask = CoBoConfigParser::agetChanBitmask(regName, channelId, itChan);
		// Compute expected bitmask
		size_t physicalChannelId = 1 + CoBoConfigParser::agetPhysicalChannelId(channelId);
		uint64_t expectedMask = 0;
		if (physicalChannelId >= 17 and physicalChannelId <= 32)
			expectedMask = threshold << (60 - 4*(physicalChannelId - 17));
		// Compare
		std::cout << "Checking that bitmask for channel " << channelId << " (physical channel " << physicalChannelId << ") in register " << regName
					<< " is " << std::hex << std::showbase << expectedMask << std::endl;
		BOOST_CHECK_EQUAL(actualMask, expectedMask);
	}
}

BOOST_AUTO_TEST_CASE(Reg8Msb_test)
{
	const std::string regName = "reg8msb";
	const uint64_t threshold = 0xF;
	// Check bitmask for individual channels
	for (size_t channelId=0; channelId < 68; channelId++)
	{
		// Create configuration
		CCfg::Object aget("Aget");
		CCfg::CConfig agetConfig(&aget);
		agetConfig("channel", channelId)("LSBThresholdValue").setIntValue(threshold);
		// Compute actual bitmask
		CCfg::CConfig::Iterator itChan = agetConfig.iterate("channel");
		itChan.next();
		uint64_t actualMask = CoBoConfigParser::agetChanBitmask(regName, channelId, itChan);
		// Compute expected bitmask
		size_t physicalChannelId = 1 + CoBoConfigParser::agetPhysicalChannelId(channelId);
		uint64_t expectedMask = 0;
		if (physicalChannelId >= 1 and physicalChannelId <= 16)
			expectedMask = threshold << (60 - 4*(physicalChannelId - 1));
		// Compare
		std::cout << "Checking that bitmask for channel " << channelId << " (physical channel " << physicalChannelId << ") in register " << regName
					<< " is " << std::hex << std::showbase << expectedMask << std::endl;
		BOOST_CHECK_EQUAL(actualMask, expectedMask);
	}
}

BOOST_AUTO_TEST_CASE(Reg9Lsb_test)
{
	const std::string regName = "reg9";
	const uint64_t threshold = 0xF;
	// Check bitmask for individual channels
	for (size_t channelId=0; channelId < 68; channelId++)
	{
		// Create configuration
		CCfg::Object aget("Aget");
		CCfg::CConfig agetConfig(&aget);
		agetConfig("channel", channelId)("LSBThresholdValue").setIntValue(threshold);
		// Compute actual bitmask
		CCfg::CConfig::Iterator itChan = agetConfig.iterate("channel");
		itChan.next();
		uint64_t actualMask = CoBoConfigParser::agetChanBitmask(regName, channelId, itChan);
		// Compute expected bitmask
		size_t physicalChannelId = 1 + CoBoConfigParser::agetPhysicalChannelId(channelId);
		uint64_t expectedMask = 0;
		if (physicalChannelId >= 33 and physicalChannelId <= 48)
			expectedMask = threshold << 4*(physicalChannelId - 33);
		// Compare
		std::cout << "Checking that bitmask for channel " << channelId << " (physical channel " << physicalChannelId << ") in register " << regName
					<< " is " << std::hex << std::showbase << expectedMask << std::endl;
		BOOST_CHECK_EQUAL(actualMask, expectedMask);
	}
}

BOOST_AUTO_TEST_CASE(Reg9Msb_test)
{
	const std::string regName = "reg9msb";
	const uint64_t threshold = 0xF;
	// Check bitmask for individual channels
	for (size_t channelId=0; channelId < 68; channelId++)
	{
		// Create configuration
		CCfg::Object aget("Aget");
		CCfg::CConfig agetConfig(&aget);
		agetConfig("channel", channelId)("LSBThresholdValue").setIntValue(threshold);
		// Compute actual bitmask
		CCfg::CConfig::Iterator itChan = agetConfig.iterate("channel");
		itChan.next();
		uint64_t actualMask = CoBoConfigParser::agetChanBitmask(regName, channelId, itChan);
		// Compute expected bitmask
		size_t physicalChannelId = 1 + CoBoConfigParser::agetPhysicalChannelId(channelId);
		uint64_t expectedMask = 0;
		if (physicalChannelId >= 49 and physicalChannelId <= 64)
			expectedMask = threshold << 4*(physicalChannelId - 49);
		// Compare
		std::cout << "Checking that bitmask for channel " << channelId << " (physical channel " << physicalChannelId << ") in register " << regName
					<< " is " << std::hex << std::showbase << expectedMask << std::endl;
		BOOST_CHECK_EQUAL(actualMask, expectedMask);
	}
}

} /* namespace rc */
} /* namespace get */
