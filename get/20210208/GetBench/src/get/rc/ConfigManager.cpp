/**
 * @file ConfigManager.cpp
 * @date 9 oct. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetBench software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "ConfigManager.h"
#include "ConfigFileManager.h"
#include "SubConfigIdSet.h"
#include "ConfigIdSet.h"

#include "CCfg/Io/BufferStorage.h"
#include "CCfg/Io/XmlFormat.h"
#include "CCfg/CConfig.h"
#include <CCfg/Object.h>
#include "CCfg/View/InheritersIterator.h"
#include "utl/Logging.h"
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/case_conv.hpp>
namespace ba = boost::algorithm;
#include <iterator>
#include <algorithm>
#include <sstream>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
ConfigManager::ConfigManager()
{
}
//_________________________________________________________________________________________________
/**
 * Returns list of triplets identifying available configuration sets.
 */
void ConfigManager::getConfigIds(ConfigIdSet & configIds) const
{
	// Get 3 sets of configuration sub-IDs
	SubConfigIdSet descrIds, prepIds, confIds;
	getSubConfigIds(ConfigId::Describe, descrIds);
	//getSubConfigIds(ConfigId::Prepare, prepIds);
	getSubConfigIds(ConfigId::Configure, confIds);

	// Constructs all triplets
	configIds.clear();
	/*SubConfigIdSet::const_iterator iter1, iter2, iter3;
	for (iter1 = descrIds.begin(); iter1 != descrIds.end(); iter1++)
	{
		for (iter2 = prepIds.begin(); iter2 != prepIds.end(); iter2++)
		{
			for (iter3 = confIds.begin(); iter3 != confIds.end(); iter3++)
			{
				configIds.insert(*iter1, *iter2, *iter3);
			}
		}
	}*/
	SubConfigIdSet::const_iterator iter1, iter2;
	for (iter1 = descrIds.begin(); iter1 != descrIds.end(); iter1++)
	{
		for (iter2 = confIds.begin(); iter2 != confIds.end(); iter2++)
		{
			configIds.insert(*iter1, *iter2, *iter2);
		}
	}
}
//_________________________________________________________________________________________________
std::auto_ptr< CCfg::Object > ConfigManager::getDescribeConfig(const std::string & subConfigId) const
{
	return getSubConfig(ConfigId::Describe, subConfigId);
}
//_________________________________________________________________________________________________
std::auto_ptr< CCfg::Object > ConfigManager::getPrepareConfig(const std::string & subConfigId) const
{
	//return getSubConfig(ConfigId::Prepare, subConfigId);
	return getSubConfig(ConfigId::Configure, subConfigId);
}
//_________________________________________________________________________________________________
std::auto_ptr< CCfg::Object > ConfigManager::getConfigureConfig(const std::string & subConfigId) const
{
	return getSubConfig(ConfigId::Configure, subConfigId);
}
//_________________________________________________________________________________________________
std::auto_ptr< CCfg::Object > ConfigManager::getHardwareDescription(const std::string & subConfigId) const
{
	return getSubConfig(ConfigId::HardwareDescr, subConfigId);
}
//_________________________________________________________________________________________________
/**
 * Returns complete configuration ID (resolving inner references).
 */
ConfigId ConfigManager::getFullConfigId(const ConfigId & id) const
{
	// Copy sub-config IDs
	ConfigId fullId(id);

	// Add inner references to hardware descriptions
	ConfigId::SubConfigIdMap::const_iterator subIdIter;
	for (subIdIter = id.subConfigIds().begin(); subIdIter != id.subConfigIds().end(); ++subIdIter)
	{
		// Loop over master description(s)
		if (subIdIter->first != ConfigId::Describe) continue;
		// Get contents
		std::auto_ptr< CCfg::Object > descrObj = getDescribeConfig(subIdIter->second);
		CCfg::CConfig config(descrObj.get());
		// Iterate
		CCfg::View::InheritersIterator inheritersIter = config("Node", "*")("Instance", "*").getInheriterAttributes();
		while (inheritersIter.hasNext())
		{
			// do not forget to call next before accessing the first element.
			CCfg::Object & nodeObj = inheritersIter.next()->toObject();
			if (nodeObj.isDefaultConfig())
					continue;
			// Get ID of hardware description
			const std::string hardwareId = CCfg::CConfig(&nodeObj)("HardwareConfigId").getStringValue();
			fullId.insert(ConfigId::HardwareDescr, hardwareId);
		}
	}

	return fullId;
}
//_________________________________________________________________________________________________
std::string ConfigManager::addDescribeConfig(std::auto_ptr< CCfg::Object > cfg, const std::string & subConfigId) const
{
	return addSubConfig(ConfigId::Describe, cfg, subConfigId);
}
//_________________________________________________________________________________________________
std::string ConfigManager::addSubConfig(const ConfigId::SubType subConfigType, const std::string & subConfig, const std::string & subConfigId) const
{
	// Create configuration from string
	CCfg::Io::BufferStorage buffer;
	buffer.write(subConfig);
	CCfg::Io::XmlFormat xml;
	std::auto_ptr<CCfg::Object> cfg = xml.parse(buffer, CCfg::Io::NoOpIncludeHandler::getInstance());

	return addSubConfig(subConfigType, cfg, subConfigId);
}
//_________________________________________________________________________________________________
/**
 * Retrieves a configuration and returns its XML representation.
 * @param subConfigType Configuration type.
 * @param subConfigId Configuration identifier.
 * @param[out] xmlConfig String containing the XML representation of the configuration.
 */
void ConfigManager::getSubConfig(const ConfigId::SubType subConfigType, const std::string & subConfigId, std::string & xmlConfig) const
{
	ConfigId::SubType actualSubType = (ConfigId::Prepare == subConfigType)?ConfigId::Configure:subConfigType;
	std::auto_ptr< CCfg::Object > cfg = getSubConfig(actualSubType, subConfigId);
	// Format config to string
	CCfg::Io::BufferStorage buffer;
	CCfg::Io::XmlFormat xml(false);
	xml.format(cfg.get(), buffer, CCfg::Io::NoOpIncludeHandler::getInstance());
	xmlConfig.clear();
	buffer.readAll(xmlConfig);
}
//_________________________________________________________________________________________________
/**
 * Returns an XML element of type ConfigSet containing all the sub-configurations corresponding to the given ID.
 * @param configId Configuration identifier.
 */
std::string ConfigManager::getConfigSet(const ConfigId & configId) const
{
	ConfigId fullId = getFullConfigId(configId);
	std::ostringstream oss;
	oss << "<ConfigSet>";
	// Loop over sub-config IDs
	ConfigId::SubConfigIdMap::const_iterator subIdIter;
	for (subIdIter = fullId.subConfigIds().begin(); subIdIter != fullId.subConfigIds().end(); ++subIdIter)
	{
		const ConfigId::SubType subType = subIdIter->first;
		const std::string subId = subIdIter->second;
		const std::string filename = ConfigFileManager::getPrefix(subIdIter->first) + subId + ".xcfg";
		oss << "\n<File filename=\"" << filename << "\" format=\"xcfg\" subConfigType=\"" << ConfigId::subTypeToString(subType) << "\" subConfigId=\"" << subId << "\">\n";
		std::string contents;
		getSubConfig(subType, subId, contents);
		oss << contents;
		oss << "</File>";
	}
	oss << "\n</ConfigSet>";

	return oss.str();
}
//_________________________________________________________________________________________________
std::string ConfigManager::addPrepareConfig(std::auto_ptr< CCfg::Object > cfg, const std::string & subConfigId) const
{
	return addSubConfig(ConfigId::Prepare, cfg, subConfigId);
}
//_________________________________________________________________________________________________
std::string ConfigManager::addConfigureConfig(std::auto_ptr< CCfg::Object > cfg, const std::string & subConfigId) const
{
	return addSubConfig(ConfigId::Configure, cfg, subConfigId);
}
//_________________________________________________________________________________________________
std::string ConfigManager::addHardwareConfig(std::auto_ptr< CCfg::Object > cfg, const std::string & subConfigId) const
{
	return addSubConfig(ConfigId::HardwareDescr, cfg, subConfigId);
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
