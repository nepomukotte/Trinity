#include "SiteItem.h"
#include "ConnectionItem.h"
#include "SchemaItem.h"
#include "RunItem.h"
#include "CCfg/Io/DatabaseStorage.h"
#include "CCfg/Io/Db/ConfigDao.h"
#include "CCfg/Io/Db/Run.h"

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::SiteItem class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
SiteItem::SiteItem(SchemaItem* parent, const Io::Db::SitePtr& site)
	: parent(parent), site(site), runs(), upToDate(false)
{
}

/**
 * Destructor.
 */
SiteItem::~SiteItem()
{
}

/**
 * Creates a new run in this site.
 */
void SiteItem::createRun()
{
	Io::Db::ConfigDao& dao = getSchema()->getDao();
	if(not dao.checkHasAllTables())
	{
		dao.replaceAllTables();
	}
	Io::Db::Run run;
	run.id = dao.allocIdRangeForRuns(1) + 1;
	run.fkSite = getEntity().getId();
	dao.insert(run);
}

/*
 * (see inherited doc)
 */
bool SiteItem::canSetName() const
{
	return true;
}

/*
 * (see inherited doc)
 */
void SiteItem::setName(const std::string& newName)
{
	parent->getDao().rename(*site, newName);
}


/*
 * (see inherited doc)
 */
const std::string& SiteItem::getName() const
{
	return site->name;
}

/*
 * (see inherited doc)
 */
size_t SiteItem::getItemCount()
{
	ensureUpToDate();
	return runs.size();
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* SiteItem::getItem(size_t index)
{
	return &runs[index];
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* SiteItem::getParent()
{
	return parent;
}

/*
 * (see inherited doc)
 */
void SiteItem::refresh()
{
	upToDate = false;
}

/*
 * (see inherited doc)
 */
void SiteItem::removeItem(size_t index, RemoveOptions options)
{
	if(options & FromDataStore)
	{
		Io::Db::ConfigDao& dao = parent->getDao();
		dao.remove(runs[index].getEntity());
	}
	if(options & FromList)
	{
		runs.erase(runs.begin() + index);
	}
}

/*
 * (see inherited doc)
 */
void SiteItem::clear()
{
	runs.clear();
}

/*
 * (see inherited doc)
 */
std::string SiteItem::buildUrl() const
{
	ConnectionItem* connec = parent->getConnection();
	return Io::DatabaseStorage::formatUrl(connec->getDbmsType(),
			connec->getUser(),
			std::string(),
			connec->getHost(),
			connec->getPort(),
			connec->getSID(),
			parent->getName(), // schema name
			site->name);
}

/*
 * (see inherited doc)
 */
SiteItem::ItemType SiteItem::getType() const
{
	return Site;
}

void SiteItem::ensureUpToDate()
{
	using namespace Io::Db;
	if(not upToDate)
	{
		site->getDAO()->setWorkingDatabase(parent->getName());
		RunPtrSetPtr siteRuns = site->getRuns();
		RunIdIdx& runIds = siteRuns->get<CCfg::Io::Db::Id>();
		for(RunIdIdx::iterator curRun = runIds.begin(),
			endRun = runIds.end(); curRun != endRun; ++curRun)
		{
			runs.push_back(new RunItem(this, *curRun));
		}
		upToDate = true;
	}
}

/**
 * Returns the DAO.
 */
Io::Db::ConfigDao& SiteItem::getDao()
{
	return parent->getDao();
}

}
}
