// =============================================================================
// Config.cpp
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "Config.h"
#include "CCfg/Scalar.h"
#include "CCfg/Attribute.h"
#include "CCfg/Enum.h"
#include "CCfg/Interval.h"

namespace CCfg
{

IoTestConfig::IoTestConfig() : daq("Daq"), config("Conf"), ambiguityTest("AmbiguityTest")
{
	Attribute& frameDuration = daq.appendScalar("frameDuration");
	frameDuration = 100;
	frameDuration.setUnit("ms");
	daq("frameDuration").setPreComment(Comment("The data acquistion duration for a frame"));
	daq("frameDuration").setPostComment(Comment("10 per second"));
	Attribute& readOutFreq = daq.appendScalar("readOutFreq");
	readOutFreq = 25;
	readOutFreq.setUnit("Mb/s");
	readOutFreq.setRange(Interval(10, 80));

	daq.appendObject("Channel", "*");
	daq.appendObject("Channel", "0");
	daq.appendObject("Channel", "1");
	daq.appendObject("Channel", "2");

	Attribute& channel = daq("Channel", "*");
	channel.appendScalar("isActive") = true;
	channel.appendObject("Flow", "SPE");
	channel("Flow", "SPE").setPreComment(Comment("Data flow parameters for each type of data"));
	channel("Flow", "SPE").setPostComment(Comment("Single Photo Electron"));

	Attribute& eventRateSPE = channel("Flow", "SPE").appendScalar("eventRate");
	eventRateSPE = 100;
	eventRateSPE.setUnit("kHz");
	Enum eventRateEnum;
	eventRateEnum = 2,100;
	eventRateSPE.setRange(eventRateEnum);

	Attribute& eventSizeSPE = channel("Flow", "SPE").appendScalar("eventSize");
	eventSizeSPE = 6;
	eventSizeSPE.setUnit("Bytes");
	Enum eventSizeEnum;
	eventSizeEnum = 6,250;
	eventSizeSPE.setRange(eventSizeEnum);

	channel.appendObject("Flow", "WFA");
	Attribute& eventRateWFA = channel("Flow", "WFA").appendScalar("eventRate");
	eventRateWFA = 2;
	Attribute& eventSizeWFA = channel("Flow", "WFA").appendScalar("eventSize");
	eventSizeWFA = 250;

	daq("Channel", "1").appendScalar("isActive") = false;


	// initialize the second configuration
	Attribute& attr = config.appendScalar("scalar");
	attr.setUnit("Tabs\tand even newlines\nin this unit!!");
	attr.setStringValue("Entities: <>&\'\". Other characters:\x07\x01\x0F");

	// initialize the third configuration
	Attribute& intString = ambiguityTest.appendScalar("intString");
	intString = "123";
	Attribute& hexString = ambiguityTest.appendScalar("hexString");
	hexString = "0xABC";
	Attribute& floatString = ambiguityTest.appendScalar("floatString");
	floatString = "5.0";
	Attribute& boolString = ambiguityTest.appendScalar("boolString");
	boolString = "true";
	ambiguityTest.appendScalar("emptyString");
	ambiguityTest.appendObject("emptyObject");
	Attribute& intReal = ambiguityTest.appendScalar("intReal");
	intReal = 5.0;
}

IoTestConfig::~IoTestConfig()
{

}

}
