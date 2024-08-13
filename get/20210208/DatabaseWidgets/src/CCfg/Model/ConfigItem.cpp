#include "ConfigItem.h"
#include "ConnectionItem.h"
#include "SchemaItem.h"
#include "VersionItem.h"
#include "CCfg/Io/DatabaseStorage.h"
#include "CCfg/Io/Db/ConfigDao.h"

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigItem class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConfigItem::ConfigItem(SchemaItem* parent, const Io::Db::ConfigPtr& cfg)
	: parent(parent), config(cfg), versions(), upToDate(false)
{
}

/**
 * Destructor.
 */
ConfigItem::~ConfigItem()
{
}

/*
 * (see inherited doc)
 */
bool ConfigItem::canSetName() const
{
	return true;
}

/*
 * (see inherited doc)
 */
void ConfigItem::setName(const std::string& newName)
{
	parent->getDao().rename(*config, newName);
}


/*
 * (see inherited doc)
 */
const std::string& ConfigItem::getName() const
{
	return config->name;
}

/*
 * (see inherited doc)
 */
size_t ConfigItem::getItemCount()
{
	ensureUpToDate();
	return versions.size();
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* ConfigItem::getItem(size_t index)
{
	return &versions[index];
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* ConfigItem::getParent()
{
	return parent;
}

/*
 * (see inherited doc)
 */
void ConfigItem::refresh()
{
	upToDate = false;
}

/*
 * (see inherited doc)
 */
void ConfigItem::removeItem(size_t index, RemoveOptions options)
{
	if(options & FromDataStore)
	{
		Io::Db::ConfigDao& dao = parent->getDao();
		dao.remove(versions[index].getEntity());
	}
	if(options & FromList)
	{
		versions.erase(versions.begin() + index);
	}
}

/*
 * (see inherited doc)
 */
void ConfigItem::clear()
{
	versions.clear();
}

/*
 * (see inherited doc)
 */
std::string ConfigItem::buildUrl() const
{
	ConnectionItem* connec = parent->getConnection();
	return Io::DatabaseStorage::formatUrl(connec->getDbmsType(),
			connec->getUser(),
			std::string(),
			connec->getHost(),
			connec->getPort(),
			connec->getSID(),
			parent->getName(), // schema name
			config->name);
}

/*
 * (see inherited doc)
 */
ConfigItem::ItemType ConfigItem::getType() const
{
	return Config;
}

void ConfigItem::ensureUpToDate()
{
	using namespace Io::Db;
	if(not upToDate)
	{
		config->getDAO()->setWorkingDatabase(parent->getName());
		ConfigVersionPtrSetPtr cfgVersions = config->getVersions();
		ConfigVersionDateIdx& versionIdx = cfgVersions->get<Date>();
		UID lastStruct = NULL_UID;
		bool structColor = true;
		for(ConfigVersionDateIdx::iterator curVersion = versionIdx.begin(),
			endVersion = versionIdx.end(); curVersion != endVersion; ++curVersion)
		{
			UID curStruct = (**curVersion).fkStructure;
			if(curStruct != lastStruct)
			{
				structColor = !structColor;
				lastStruct = curStruct;
			}
			versions.push_back(new VersionItem(this, *curVersion, structColor));
		}
		upToDate = true;
	}
}

}
}
