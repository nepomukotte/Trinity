#include "ConfigDao.h"
#include "Config.h"
#include "InsertionProcess.h"
#include "AbstractUpdateController.h"
#include "CCfg/Attribute.h"
#include <ctime>
#include <cassert>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::Config class methods body
////////////////////////////////////////////////////////////////////////////////

const char Config::TABLE_NAME[] = "CCfg_Config";

/**
 * Default constructor.
 */
Config::Config()
{
}

/**
 * Destructor.
 */
Config::~Config()
{
}

/**
 * Equality operator.
 */
bool Config::operator==(const Config& entity)
{
	return id == entity.id and
		name == entity.name;
}

/**
 * @fn bool Config::operator!=(const Config& entity)
 * Unequality operator.
 */

/**
 * Inserts a in-memory configuration as a new version of this configuration.
 * @param configRoot the root of the configuration.
 * @param controller monitors and controls structure update.
 * This object will receive notification about structural differences and will
 * have to accept or refuse the whole structure update.
 * @return Returns the newly created version, or null if the controller canceled
 * the structure update.
 */
ConfigVersionPtr Config::insertNewVersion(const Attribute* configRoot, AbstractUpdateController& controller)
{
	ConfigVersionPtr newVersion(0);
	InsertionProcess inserter(this, const_cast<Attribute*>(configRoot), controller);
	if(inserter.insert())
	{
		newVersion = inserter.getNewVersion();
	}
	return newVersion;
}

/**
 * Returns the list of all classes used for this configuration.
 */
ConfigClassPtrSetPtr Config::getClasses()
{
	assert(dao != 0);
	return dao->getClassesOfConfig(*this);
}

/**
 * Returns the list of structure versions of this configuration.
 */
ConfigStructureVersionPtrSetPtr Config::getStructureVersions()
{
	assert(dao != 0);
	return dao->getStructureVersions(*this);
}

/**
 * Returns the list of versions of this configuration.
 */
ConfigVersionPtrSetPtr Config::getVersions()
{
	assert(dao != 0);
	return dao->getVersionsOfConfig(*this);
}

/**
 * Returns the most recent structure version of this configuration.
 */
ConfigStructureVersionPtr Config::getLastStructureVersion()
{
	assert(dao != 0);
	return dao->getLastStructureVersion(*this);
}

/**
 * Returns the structure version matching exactly the specified date.
 */
ConfigStructureVersionPtr Config::getStructureVersionAt(const std::tm& date)
{
	assert(dao != 0);
	return dao->getStructureVersionAt(date);
}

/**
 * Returns the most recent version of this configuration.
 */
ConfigVersionPtr Config::getLastVersion()
{
	assert(dao != 0);
	return dao->getLastVersion(*this);
}

/**
 * Returns the version record matching exactly the specified date.
 */
ConfigVersionPtr Config::getVersionAt(const std::tm& date)
{
	assert(dao != 0);
	return dao->getVersionAt(date);
}

/**
 * Instanciates a new structure version entity for this configuration
 */
ConfigStructureVersionPtr Config::createNewStructure()
{
	assert(dao != 0);
	ConfigStructureVersionPtr sversion(new ConfigStructureVersion);
	sversion->setDAO(dao);
	sversion->id = dao->allocIdRangeForStructureVersions(1) + 1;
	sversion->fkConfig = this->id;
	time_t currentTime = time(0);
	sversion->creationDate = *gmtime(&currentTime);
	return sversion;
}

/**
 * Creates a new Config entity (it must still be inserted manually).
 */
ConfigPtr Config::createConfig(const std::string& cfgName, ConfigDao* dao)
{
	assert(dao != 0);
	ConfigPtr config(new Config);
	config->setDAO(dao);
	config->id = dao->allocIdRangeForConfigs(1) + 1;
	config->name = cfgName;
	return config;
}


}
}
}
