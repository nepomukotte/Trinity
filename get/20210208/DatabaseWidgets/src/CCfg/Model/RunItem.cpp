#include "RunItem.h"
#include "SiteItem.h"
#include "ConditionItem.h"
#include "SchemaItem.h"
#include "ConnectionItem.h"
#include "VersionItem.h"
#include "CCfg/Io/Db/RunCondition.h"
#include "Utils/Date/StdTmHelpers.h"
#include <sstream>
#include <stdexcept>

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::RunItem class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
RunItem::RunItem(SiteItem* parent, const Io::Db::RunPtr& run)
	: parent(parent), run(run), conditions(), name(), upToDate(false)
{
	updateDisplayName();
}

/**
 * Destructor.
 */
RunItem::~RunItem()
{
}

/*
 * (see inherited doc)
 */
bool RunItem::isLeaf() const
{
	return false;
}

/*
 * (see inherited doc)
 */
bool RunItem::canSetName() const
{
	return false;
}

/*
 * (see inherited doc)
 */
void RunItem::setName(const std::string& )
{
	throw std::runtime_error("cannot change the name of a run");
}

/*
 * (see inherited doc)
 */
const std::string& RunItem::getName() const
{
	return name;
}

/*
 * (see inherited doc)
 */
const std::string& RunItem::getDisplayName() const
{
	return name;
}

/*
 * (see inherited doc)
 */
size_t RunItem::getItemCount()
{
	ensureUpToDate();
	return conditions.size();
}

void RunItem::ensureUpToDate()
{
	using namespace Io::Db;
	if (not upToDate)
	{
		run->getDAO()->setWorkingDatabase(getSite()->getSchema()->getName());
		RunConditionPtrSetPtr runConditions = run->getConditions();
		RunConditionPtrIdIdx& conditionIdx = runConditions->get<CCfg::Io::Db::Id>();
		for (RunConditionPtrIdIdx::iterator curCondition = conditionIdx.begin(),
			endCondition = conditionIdx.end(); curCondition != endCondition; ++curCondition)
		{
			conditions.push_back(new ConditionItem(this, *curCondition));
		}
		upToDate = true;
	}
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* RunItem::getItem(size_t index)
{
	return &conditions[index];
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* RunItem::getParent()
{
	return parent;
}

/*
 * (see inherited doc)
 */
void RunItem::refresh()
{
	upToDate = false;
}

/*
 * (see inherited doc)
 */
void RunItem::removeItem(size_t index, RemoveOptions options)
{
	if(options & FromDataStore)
	{
		Io::Db::ConfigDao& dao = parent->getDao();
		dao.remove(conditions[index].getEntity());
	}
	if(options & FromList)
	{
		conditions.erase(conditions.begin() + index);
	}
}

/*
 * (see inherited doc)
 */
void RunItem::clear()
{
	conditions.clear();
}

/*
 * (see inherited doc)
 */
std::string RunItem::buildUrl() const
{
	ConnectionItem* connec = parent->getSchema()->getConnection();
	std::ostringstream runStr;
	runStr << parent->getName() << '#' << run->number;
	return Io::DatabaseStorage::formatUrl(connec->getDbmsType(),
			connec->getUser(),
			std::string(),
			connec->getHost(),
			connec->getPort(),
			connec->getSID(),
			parent->getSchema()->getName(),
			runStr.str());
}

/*
 * (see inherited doc)
 */
RunItem::ItemType RunItem::getType() const
{
	return Run;
}

void RunItem::updateDisplayName()
{
	std::ostringstream nameStream;
	nameStream << "Run" << run->number;
	name = nameStream.str();
}

void RunItem::createCondition(const VersionItem* versionItem)
{
	CCfg::Io::Db::RunConditionPtr condition = run->createNewCondition(versionItem->getEntity());
	parent->getSchema()->getDao().insert(*condition);
	upToDate = false;
}

}
}
