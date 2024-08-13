/**
 * @file ConfigDatabaseManager.cpp
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
#include "ConfigDatabaseManager.h"
#include <CCfg/Io/Document.h>
#include <CCfg/Io/FormatFactory.h>
#include <CCfg/Io/AbstractFormatFactory.h>
#include <CCfg/Io/DatabaseFormatFactory.h>
#include <CCfg/Io/StorageFactory.h>
#include <CCfg/Io/DatabaseStorageFactory.h>
#include <CCfg/Io/Db/ConfigVersion.h>
#include <Utils/Date/StdTmHelpers.h>
#include <CCfg/Io/Db/Config.h>
#include "utl/Logging.h"
using namespace CCfg::Io;
using namespace CCfg::Io::Db;
#include <sstream>
#include <iomanip>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
bool ConfigDatabaseManager::InvalidLoginCallback::requestLogin(CCfg::Io::DbmsType /* dbms */,
		const std::string& /* host */, const std::string& /* port */,
		const std::string& /* sid */, std::string& /* user */, std::string& /* pwd */)
{
	return false;
}
//_________________________________________________________________________________________________
bool ConfigDatabaseManager::dbFactoriesReady = false;
//_________________________________________________________________________________________________
ConfigDatabaseManager::ConfigDatabaseManager(const std::string & repoUrl) : ConfigManager()
{
	initFactories();
	// Store connection parameters
	dbUrl = Url(repoUrl);
	user = dbUrl.getUser();
	password = dbUrl.getPassword();

	LOG_INFO() << "Creating manager for configuration database '" << obfuscateUrl(dbUrl) << '\'';
}
//_________________________________________________________________________________________________
CCfg::Io::Url ConfigDatabaseManager::obfuscateUrl(const CCfg::Io::Url & url)
{
	CCfg::Io::Url obfUrl(url);
	obfUrl.setUser("");
	obfUrl.setPassword("");
	return obfUrl;
}
//_________________________________________________________________________________________________
void ConfigDatabaseManager::initFactories()
{
	if (not dbFactoriesReady)
	{
		DatabaseFormatFactory* databaseFormatFactory = new DatabaseFormatFactory;
		std::auto_ptr<AbstractFormatFactory> formatFactory(databaseFormatFactory);
		databaseFormatFactory->setUpdateController(UpdateControllerPtr(new Db::DefaultUpdateController));
		FormatFactory::instance().registerFormat(formatFactory);

		DatabaseStorageFactory* databaseStorageFactory = new DatabaseStorageFactory;
		std::auto_ptr<AbstractStorageFactory> storageFactory(databaseStorageFactory);
		databaseStorageFactory->setLoginCallback(LoginCallbackPtr(new InvalidLoginCallback));
		StorageFactory::instance().registerStorage(storageFactory);

		dbFactoriesReady = true;
	}
}
//_________________________________________________________________________________________________
std::string ConfigDatabaseManager::getConfigName(const ConfigId::SubType type) const
{
	switch (type)
	{
	case ConfigId::Describe:
		return "description";
	case ConfigId::Prepare:
		return "prepare";
	case ConfigId::Configure:
		return "configuration";
	case ConfigId::HardwareDescr:
		return "hardwareDescription";
	default:
		LOG_ERROR() << "Invalid configuration type '" << type << '\'';
		break;
	}
	return "";
}
//_________________________________________________________________________________________________
CCfg::Io::Url ConfigDatabaseManager::formatUrl(const ConfigId::SubType type, const std::string & configId) const
{
	// Build URL
	Url configUrl = dbUrl;
	configUrl.appendPath(getConfigName(type));
	configUrl.addQueryItem("date", configId);

	return configUrl;
}
//_________________________________________________________________________________________________
std::auto_ptr< CCfg::Object > ConfigDatabaseManager::getConfig(const CCfg::Io::Url & cfgUrl) const
{
	CCfg::Io::Url configUrl(cfgUrl);
	LOG_DEBUG() << "Loading configuration with URL '" << configUrl << "' ...";

	configUrl.setUser(user);
	configUrl.setPassword(password);
	CCfg::Io::Document doc;
	try
	{
		doc.loadFromUrl(configUrl);
	}
	catch (const CCfg::Io::Document::IoError & e)
	{
		LOG_ERROR() << "Could not retrieve the configuration at '" << obfuscateUrl(configUrl) << " (check your credentials?).";
		throw ConfigNotFound(configUrl.format());
	}

	return doc.detachConfig();
}
//_________________________________________________________________________________________________
std::auto_ptr< CCfg::Object > ConfigDatabaseManager::getSubConfig(const ConfigId::SubType type, const std::string & configId) const
{
	// Check configuration exists
	ConfigPtr config = getConfigDao().searchConfig(getConfigName(type));
	if (not config)
		throw ConfigIdNotFound(configId);

	// Check version exists
	std::istringstream iss(configId);
	std::tm date;
	iss >> date;
	if (not iss)
		throw ConfigIdNotFound(configId);
	ConfigVersionPtr version = config->getVersionAt(date);
	if (not version)
		throw ConfigIdNotFound(configId);

	// Build URL
	CCfg::Io::Url configUrl = formatUrl(type, configId);

	// Read configuration
	try
	{
		return getConfig(configUrl);
	}
	catch (const ConfigNotFound & e)
	{
		throw ConfigIdNotFound(configId);
	}
	return std::auto_ptr< CCfg::Object >();
}
//_________________________________________________________________________________________________
Db::ConfigDao & ConfigDatabaseManager::getConfigDao() const
{
	static DatabaseStorage* dbStorage = 0;
	static std::auto_ptr<CCfg::Io::Storage> storage;

	if (dbStorage == 0)
	{
		try
		{
			// Get DAO
			Url connectionUrl(dbUrl);
			connectionUrl.setUser(user);
			connectionUrl.setPassword(password);
			storage = StorageFactory::instance().createStorage(connectionUrl);
			if (storage.get())
			{
				dbStorage = dynamic_cast<DatabaseStorage*>(storage.get());
				dbStorage->open();
			}
		}
		catch (const std::exception & e)
		{
			LOG_ERROR() << e.what();
			throw ConnectionError(dbUrl.format(), e.what());
		}
	}

	if (dbStorage == 0 or not dbStorage->isOpened())
		throw ConnectionError(dbUrl.format());

	return dbStorage->getDao();
}
//_________________________________________________________________________________________________
void ConfigDatabaseManager::getSubConfigIds(const ConfigId::SubType type, SubConfigIdSet & configIds) const
{
	// Clear list
	configIds.clear();

	// Check contents of database
	Db::ConfigDao & dao = getConfigDao();

	try
	{
		ConfigPtr config = dao.searchConfig(getConfigName(type));
		if (config)
		{
			ConfigVersionPtrSetPtr versions = config->getVersions();
			ConfigVersionDateIdx & dateIdx = versions->get<Date>();
			for (ConfigVersionDateIdx::const_iterator versionIt = dateIdx.begin(); versionIt != dateIdx.end(); ++versionIt)
			{
				ConfigVersion & version = **versionIt;
				std::ostringstream oss;
				oss << version.getCreationDate();
				configIds.insert(oss.str());
			}
		}
	}
	catch (const std::exception & e)
	{
		throw ConnectionError(dbUrl.format(), e.what());
	}
}
//_________________________________________________________________________________________________
std::string ConfigDatabaseManager::addSubConfig(const ConfigId::SubType type, std::auto_ptr< CCfg::Object > cfg, const std::string & /* configId */) const
{
	try
	{
		// Check configuration exists
		ConfigPtr config = getConfigDao().searchConfig(getConfigName(type));
		if (not config)
		{
			config = CCfg::Io::Db::Config::createConfig(getConfigName(type), &getConfigDao());
		}

		// Create new version
		CCfg::Io::Db::DefaultUpdateController updateController;
		ConfigVersionPtr version = config->insertNewVersion(cfg.get(), updateController);

		std::ostringstream oss;
		oss << version->creationDate;
		const std::string configId = oss.str();
		LOG_INFO() << "Added version '" << configId << "' of configuration '" << config->name << '\'';
		return configId;
	}
	catch (const std::exception & e)
	{
		throw ConnectionError(dbUrl.format(), e.what());
	}
	return "";
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
