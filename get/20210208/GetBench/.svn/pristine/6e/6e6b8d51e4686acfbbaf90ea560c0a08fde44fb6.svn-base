/**
 * @file ConfigFileManager.cpp
 * @date 11 oct. 2012
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

#include "ConfigFileManager.h"
#include "ConfigIdSet.h"
#include "SubConfigIdSet.h"

#include <CCfg/Io/Document.h>
#include <CCfg/Io/FileStorage.h>
#include <CCfg/Io/XmlFormat.h>
#include "utl/Logging.h"

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <boost/algorithm/string/predicate.hpp>
#include <algorithm>
#include <iomanip>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
ConfigFileManager::ConfigFileManager(const std::string & repoUrl)
	: ConfigManager(), directoryUrl(repoUrl)
{
	if (directoryUrl == ".")
		directoryUrl = fs::current_path().string();
	LOG_DEBUG() << "Creating manager for configuration files in '" << directoryUrl << "'";
}
//_________________________________________________________________________________________________
std::auto_ptr< CCfg::Object > ConfigFileManager::getConfig(const CCfg::Io::Url & cfgUrl) const
{
	LOG_DEBUG() << "Reading configuration file '" << cfgUrl << "'";
	// Read configuration
	CCfg::Io::Document doc;
	doc.loadFromFile(cfgUrl.format(), false);
	return doc.detachConfig();
}
//_________________________________________________________________________________________________
std::string ConfigFileManager::getPrefix(const ConfigId::SubType type)
{
	switch (type)
	{
	case ConfigId::Describe:
		return "describe-";
	case ConfigId::Prepare:
		return "prepare-";
	case ConfigId::Configure:
		return "configure-";
	case ConfigId::HardwareDescr:
		return "hardwareDescription_";
	default:
		LOG_ERROR() << "Invalid configuration type '" << type << '\'';
		break;
	}
	return "";
}
//_________________________________________________________________________________________________
std::auto_ptr< CCfg::Object > ConfigFileManager::getSubConfig(const ConfigId::SubType subType, const std::string & subConfigId) const
{
	// Build path
	fs::path configPath = fs::path(directoryUrl) / fs::path(getPrefix(subType) + subConfigId);
	configPath.replace_extension(".xcfg");

	// Check file exists
	if (not fs::exists(configPath))
	{
		configPath.replace_extension(".xml");
		if (not fs::exists(configPath))
		{
			LOG_DEBUG() << "Configuration file '" << configPath.replace_extension(".xcfg").string() << "' does not exist!";
			throw ConfigIdNotFound(subConfigId);
		}
	}

	return getConfig(CCfg::Io::Url(configPath.string()));
}
//_________________________________________________________________________________________________
void ConfigFileManager::getSubConfigIds(const ConfigId::SubType type, SubConfigIdSet & subConfigIds) const
{
	// Clear list
	subConfigIds.clear();

	// Check contents of directory
	fs::path directoryPath(directoryUrl);
	if (fs::exists(directoryPath) and fs::is_directory(directoryPath))
	{
		typedef std::vector<fs::path> PathSet;
		PathSet paths;
		std::copy(fs::directory_iterator(directoryPath), fs::directory_iterator(),
				std::back_inserter(paths));
		std::sort(paths.begin(), paths.end()); // sort, since directory iteration is not ordered on some file systems

		for (PathSet::const_iterator it = paths.begin(); it != paths.end(); ++it)
		{
#if BOOST_FILESYSTEM_VERSION >= 3
			std::string fn = it->filename().string();  // extract the filename from the path
			std::string ext = it->extension().string();  // extract the filename from the path
#else
			std::string fn = it->filename();
			std::string ext = it->extension();
#endif
			if (boost::starts_with(fn, getPrefix(type)) and (ext == ".xml" or ext == ".xcfg"))
			{
#if BOOST_FILESYSTEM_VERSION >= 3
				std::string configBasename = fs::path(*it).replace_extension().filename().string();
#else
				std::string configBasename = fs::path(*it).replace_extension("").filename();
#endif
				std::string configId = configBasename.substr(getPrefix(type).size());
				subConfigIds.insert(configId);
			}
		}
	}
}
//_________________________________________________________________________________________________
std::string ConfigFileManager::addSubConfig(const ConfigId::SubType type, std::auto_ptr< CCfg::Object > cfg, const std::string & subConfigId) const
{
	// Build path
	fs::path configPath = fs::path(directoryUrl) / fs::path(getPrefix(type) + subConfigId);
	configPath.replace_extension(".xcfg");

	// Check file exists
	if (fs::exists(configPath))
	{
		LOG_WARN() << "Overwriting configuration!";
		fs::remove(configPath);
	}

	// Write config
#if BOOST_FILESYSTEM_VERSION >= 3
	CCfg::Io::FileStorage outFile(configPath.c_str(), CCfg::Io::FileStorage::Write | CCfg::Io::FileStorage::Truncate);
#else
        CCfg::Io::FileStorage outFile(configPath.string().c_str(), CCfg::Io::FileStorage::Write | CCfg::Io::FileStorage::Truncate);
#endif
	CCfg::Io::XmlFormat xml;
	xml.format(cfg.get(), outFile, CCfg::Io::NoOpIncludeHandler::getInstance());
	outFile.close();

	return subConfigId;
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
