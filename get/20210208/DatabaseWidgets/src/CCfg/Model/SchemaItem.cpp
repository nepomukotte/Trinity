#include "SchemaItem.h"
#include "ConfigItem.h"
#include "SiteItem.h"
#include "ConnectionItem.h"
#include "CCfg/Io/Db/Config.h"
#include <stdexcept>

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::SchemaItem class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
SchemaItem::SchemaItem(ConnectionItem* parent, const Io::Db::ConfigDatabasePtr& db, DisplayMode mode)
	: parent(parent), database(db), configurations(), sites(), upToDate(false), displayMode(mode)
{
}

/**
 * Destructor.
 */
SchemaItem::~SchemaItem()
{
}

/**
 * Changes the password of the user associated with this schema.
 */
void SchemaItem::changeUserPassword(const std::string& newPassword)
{
	getDao().changeUserPassword(database->name, newPassword);
}

/**
 * Creates a new configuration in this schema.
 */
void SchemaItem::createConfig(const std::string& name)
{
	Io::Db::ConfigDao& dao = getDao();
	dao.setWorkingDatabase(database->name);
	if(not dao.checkHasAllTables())
	{
		dao.replaceAllTables();
	}
	Io::Db::Config config;
	config.id = dao.allocIdRangeForConfigs(1) + 1;
	config.name = name;
	dao.insert(config);
}

/**
 * Creates a new site in this schema.
 */
void SchemaItem::createSite(const std::string& name)
{
	Io::Db::ConfigDao& dao = getDao();
	dao.setWorkingDatabase(database->name);
	if(not dao.checkHasAllTables())
	{
		dao.replaceAllTables();
	}
	Io::Db::Site site;
	site.id = dao.allocIdRangeForSites(1) + 1;
	site.name = name;
	dao.insert(site);
}

/*
 * (see inherited doc)
 */
bool SchemaItem::canSetName() const
{
	return false;
}

/*
 * (see inherited doc)
 */
void SchemaItem::setName(const std::string& /*newName*/)
{
	throw std::runtime_error("cannot change the name of a schema");
}

/*
 * (see inherited doc)
 */
const std::string& SchemaItem::getName() const
{
	return database->name;
}

/*
 * (see inherited doc)
 */
size_t SchemaItem::getItemCount()
{
	ensureUpToDate();
	return configurations.size()+sites.size();
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* SchemaItem::getItem(size_t index)
{
	const size_t siteCount = sites.size();
	if (index < siteCount)
		return &sites[index];
	return &configurations[index-siteCount];
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* SchemaItem::getParent()
{
	return parent;
}

/*
 * (see inherited doc)
 */
void SchemaItem::refresh()
{
	upToDate = false;
}

/*
 * (see inherited doc)
 */
void SchemaItem::removeItem(size_t index, RemoveOptions options)
{
	const size_t siteCount = sites.size();

	if(options & FromDataStore)
	{
		Io::Db::ConfigDao& dao = getDao();
		if (index < siteCount)
			dao.remove(sites[index].getEntity());
		else
			dao.remove(configurations[index-siteCount].getEntity());
	}
	if(options & FromList)
	{
		if (index < siteCount)
			sites.erase(sites.begin() + index);
		else
			configurations.erase(configurations.begin() + index - siteCount);
	}
}

/*
 * (see inherited doc)
 */
void SchemaItem::clear()
{
	configurations.clear();
	sites.clear();
}

/*
 * (see inherited doc)
 */
std::string SchemaItem::buildUrl() const
{
	return Io::DatabaseStorage::formatUrl(parent->getDbmsType(),
			parent->getUser(),
			std::string(),
			parent->getHost(),
			parent->getPort(),
			parent->getSID(),
			database->name,
			std::string());
}

/*
 * (see inherited doc)
 */
SchemaItem::ItemType SchemaItem::getType() const
{
	return Schema;
}

/**
 * Returns the DAO for the parent connection, and set this schema as the current
 * one for queries.
 */
Io::Db::ConfigDao& SchemaItem::getDao()
{
	Io::Db::ConfigDao& dao = *database->getDAO();
	dao.setWorkingDatabase(database->name);
	return dao;
}

void SchemaItem::ensureUpToDate()
{
	using namespace Io::Db;
	if(not upToDate)
	{
		ConfigDao& dao = getDao();
		if(dao.checkHasAllTables())
		{
			if (displayMode & Config)
			{
				ConfigPtrSetPtr configs = dao.getConfigs();
				ConfigIdIdx& configIdx = configs->get<Io::Db::Id>();
				for(ConfigIdIdx::iterator curConfig = configIdx.begin(),
					endConfig = configIdx.end(); curConfig != endConfig; ++curConfig)
				{
					configurations.push_back(new ConfigItem(this, *curConfig));
				}
			}
			if (displayMode & Run)
			{
				SitePtrSetPtr siteSet = dao.getSites();
				SiteIdIdx& siteIdx = siteSet->get<Io::Db::Id>();
				for(SiteIdIdx::iterator curSite = siteIdx.begin(),
					endSite = siteIdx.end(); curSite != endSite; ++curSite)
				{
					sites.push_back(new SiteItem(this, *curSite));
				}
			}
		}
		upToDate = true;
	}
}

}
}
