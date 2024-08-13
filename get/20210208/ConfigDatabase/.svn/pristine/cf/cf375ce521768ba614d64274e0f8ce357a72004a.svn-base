#include "ConfigDao.h"
#include "ConfigDatabase.h"
#include <cassert>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigDatabase class methods body
////////////////////////////////////////////////////////////////////////////////


/**
 * Constructor.
 */
ConfigDatabase::ConfigDatabase()
{
}

/**
 * Destructor.
 */
ConfigDatabase::~ConfigDatabase()
{
}

/**
 * Returns the name of the database.
 */
const std::string& ConfigDatabase::getName() const
{
	return name;
}

/**
 * Returns the list of configurations in this database.
 */
ConfigPtrSetPtr ConfigDatabase::getConfigs()
{
	assert(dao != 0);
	dao->setWorkingDatabase(name);
	return dao->getConfigs();
}

/**
 * Returns the configuration bearing the specified name, if it exists in this database.
 */
ConfigPtr ConfigDatabase::searchConfig(const std::string& cfgName)
{
	assert(dao != 0);
	dao->setWorkingDatabase(name);
	return dao->searchConfig(cfgName);
}

/**
 * Returns whether the current database has all the required tables for the
 * configuration system.
 */
bool ConfigDatabase::checkHasAllTables()
{
	assert(dao != 0);
	dao->setWorkingDatabase(name);
	return dao->checkHasAllTables();
}

}
}
}
