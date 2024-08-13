/**
 * @file NodeAtTop.cpp
 * Implementation of the Class NodeAtTop, which translates
 * from the CompoundConfig structure to the MDaq structure
 * It expects a hierarchical tree Node->Device->Register->BitField,
 * in which Node is at the top of the CompoundConfig configuration file
 * Example:
 * -# Node[CoBo]
 * -# Device[AsAd0]
 *  - registerAccess
 *  - registerWidth
 *  - baseAddress
 * -# Register[led]
 *  - offset
 * -# BitField[l0]
 *  - offset
 *  - width
 *  - readOnly
 *
 */

#include "NodeAtTop.h"
#include <stdexcept>
#include <unistd.h>

using namespace mdaq::hw;
using std::runtime_error;

namespace mdaq
{
namespace cfg
{

NodeAtTop::NodeAtTop()
{
}

NodeAtTop::~NodeAtTop()
{
}

/**
 * @brief translates the node
 *
 * The first name of the CompoundConfig description should be "Setup[NameOfNode]",
 *
 * @param cfg the CompoundConfig configuration, read from a file/database beforehand.
 */
void NodeAtTop::translateNode(CCfg::CConfig& cfg)
{
	if (cfg.getIndex() == "")
	{
		throw runtime_error("Index of the MDaq node description cannot be empty.");
	}

	this->id = cfg.getIndex();
}

/**
 * @brief translates the device
 *
 * It looks for something like "Device[NameOfDevice]"
 * and translates the different "baseAddress", "registerAccess", "registerWidth"
 * into the MDaq description
 *
 * @param cfg the CompoundConfig configuration
 */
DeviceConfig& NodeAtTop::translateDevice(CCfg::CConfig& cfg)
{
	int n = this->devices.size();
	if (cfg.getIndex() == "")
	{
		throw runtime_error("Index of Device can't be empty");
	}

	if (cfg.getIndex() != "*")
	{
		this->devices.resize(n + 1);
		this->devices[n].descr.name = cfg.getIndex();
		this->devices[n].descr.baseAddress = cfg("baseAddress").getHexValue();
		this->devices[n].descr.registerAccess = cfg("registerAccess").getStringValue();
		this->devices[n].descr.registerWidth = cfg("registerWidth").getIntValue();

		return this->devices[n];
	}
	else
	{
		return deviceConfigNull;
	}

}

/**
 * @brief translates the register
 *
 * It looks for something like "Device[NameOfDevice]"
 * and translates "name", "readOnly", "offset"
 * into the MDaq description
 *
 * @param cfg The CompoundConfig configuration
 * @param dev The DeviceConfig&, which is the structure of the parent device inside MDaq
 */
RegisterConfig& NodeAtTop::translateRegister(CCfg::CConfig& cfg, DeviceConfig& dev)
{
	int n = dev.registers.size();

	if (cfg.getIndex() == "")
	{
		throw runtime_error("Index of Register can't be empty");
	}

	if (cfg.getIndex() != "*")
	{
		dev.registers.resize(n + 1);
		dev.registers[n].descr.name = cfg.getIndex();
		dev.registers[n].descr.offset = cfg("offset").getHexValue();

		return dev.registers[n];
	}
	else
	{
		return registerConfigNull;
	}
}

/**
 * @brief translates the BitField
 *
 * It looks for something like "BitField[NameOfBitField]"
 * and translates "name", "offset", "width"
 * into the MDaq description of the BitField
 *
 * @param cfg The CompoundConfig configuration
 * @param reg The RegisterConfig&, which is the structure of the parent register inside MDaq
 */
void NodeAtTop::translateBitField(CCfg::CConfig& cfg, RegisterConfig& reg)
{
	int n = reg.fields.size();

	if (cfg.getIndex() == "")
	{
		throw runtime_error("Index of BitField can't be empty");
	}

	if (cfg.getIndex() != "*")
	{
		reg.fields.resize(n + 1);
		reg.fields[n].name = cfg.getIndex();
		reg.fields[n].offset = cfg("offset").getIntValue();
		reg.fields[n].width = cfg("width").getIntValue();
	}
}

/**
 * @brief translates from the CompoundConfig structure to the MDaq structure
 *
 * It expects a hierarchical tree Node->Device->Register->BitField.
 * Example:
 * -# Node[CoBo]
 * -# Device[AsAd0]
 *  - registerAccess
 *  - registerWidth
 *  - baseAddress
 * -# Register[led]
 *  - offset
 * -# BitField[l0]
 *  - offset
 *  - width
 *  - readOnly
 *
 * @param cfg The CompoundConfig configuration
 */
void NodeAtTop::fromCCfg(const CCfg::CConfig& nodeCfg)
{
	CCfg::CConfig cfg(nodeCfg);
	translateNode(cfg);
	CCfg::CConfig::Iterator iDev = cfg.iterate("Device");
	while (iDev.hasNext())
	{
		CCfg::CConfig cfgDev = iDev.next();
		DeviceConfig& dev = translateDevice(cfgDev);
		if (dev != deviceConfigNull)
		{
			CCfg::CConfig::Iterator iReg = cfgDev.iterate("Register");
			while (iReg.hasNext())
			{
				CCfg::CConfig cfgReg = iReg.next();
				RegisterConfig& reg = translateRegister(cfgReg, dev);
				if (reg != registerConfigNull)
				{
					CCfg::CConfig::Iterator iField = cfgReg.iterate("BitField");
					while (iField.hasNext())
					{
						CCfg::CConfig cfgField = iField.next();
						translateBitField(cfgField, reg);
					}
				}
			}
		}
	}
	if (this->devices.size() == 0)
	{
		throw runtime_error("There should be at least one device");
	}
}

} //cfg
} //mdaq
