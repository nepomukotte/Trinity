/*
 * NodeAtTop.h
 *
 *  Created on: Mar 17, 2010
 *      Author: jchavas
 */

/**
* @file NodeAtTop.h
* Class NodeAtTop translating a CompoundConfig configuration into a MDaq structure.
*
*/

#ifndef mdaq_cfg_NodeAtTop_INCLUDED
#define mdaq_cfg_NodeAtTop_INCLUDED

#include "mdaq/hw/Control.h"
#include ICE_HEADER(Communicator.h)
#include ICE_HEADER(Initialize.h)
#include ICE_HEADER(Properties.h)
//#include ICE_HEADER(Exception)
#include "CCfg/CConfig.h"
#include "CCfg/Io/Document.h"

namespace mdaq
{
namespace cfg
{

/**
 * @brief Translates a CompoundConfig description into a MDaq description
 * with a Node being at the top of the structure
 *
 */
class NodeAtTop : public mdaq::hw::NodeConfig
{
public:
	NodeAtTop();
	virtual ~NodeAtTop();
	mdaq::hw::DeviceConfig deviceConfigNull;
	mdaq::hw::RegisterConfig registerConfigNull;
	void translateNode(CCfg::CConfig & cfg);
	mdaq::hw::DeviceConfig& translateDevice(CCfg::CConfig & cfg);
	mdaq::hw::RegisterConfig& translateRegister(CCfg::CConfig & cfg, mdaq::hw::DeviceConfig & dev);
	void translateBitField(CCfg::CConfig & cfg, mdaq::hw::RegisterConfig & reg);
	void fromCCfg(const CCfg::CConfig & cfg);
};

} //cfg
} //mdaq

#endif /* mdaq_cfg_NodeAtTop_INCLUDED */
