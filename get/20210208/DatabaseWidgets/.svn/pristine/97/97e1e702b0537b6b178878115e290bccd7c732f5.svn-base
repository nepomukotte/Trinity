#include "ConditionItem.h"
#include "RunItem.h"
#include "SchemaItem.h"
#include "ConnectionItem.h"
#include "SiteItem.h"
#include "Utils/Date/StdTmHelpers.h"
#include <sstream>
#include <stdexcept>

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::RunConditionItem class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConditionItem::ConditionItem(RunItem* parent, const Io::Db::RunConditionPtr& condition)
	: parent(parent), runCondition(condition), name()
{
	updateDisplayName();
}

/**
 * Destructor.
 */
ConditionItem::~ConditionItem()
{
}

/*
 * (see inherited doc)
 */
bool ConditionItem::isLeaf() const
{
	return true;
}

/*
 * (see inherited doc)
 */
bool ConditionItem::canSetName() const
{
	return false;
}

/*
 * (see inherited doc)
 */
void ConditionItem::setName(const std::string& )
{
	throw std::runtime_error("cannot change the name of a run condition");
}

/*
 * (see inherited doc)
 */
const std::string& ConditionItem::getName() const
{
	return name;
}

/*
 * (see inherited doc)
 */
const std::string& ConditionItem::getDisplayName() const
{
	return name;
}

/*
 * (see inherited doc)
 */
size_t ConditionItem::getItemCount()
{
	return 0;
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* ConditionItem::getItem(size_t /*index*/)
{
	return 0;
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* ConditionItem::getParent()
{
	return parent;
}

/*
 * (see inherited doc)
 */
void ConditionItem::refresh()
{
	updateDisplayName();
}

/*
 * (see inherited doc)
 */
void ConditionItem::removeItem(size_t /*index*/, RemoveOptions /*options*/)
{
}

/*
 * (see inherited doc)
 */
void ConditionItem::clear()
{
}

/*
 * Builds URL of configuration version.
 * (see inherited doc)
 */
std::string ConditionItem::buildUrl() const
{
	SchemaItem* schema = getRun()->getSite()->getSchema();
	ConnectionItem* connec = schema->getConnection();
	Io::Db::ConfigVersionPtr version = runCondition->getConfigVersion();
	return Io::DatabaseStorage::formatUrl(connec->getDbmsType(),
			connec->getUser(),
			std::string(),
			connec->getHost(),
			connec->getPort(),
			connec->getSID(),
			schema->getName(),
			version->getStructure()->getConfig()->name,
			&version->creationDate);
}

/*
 * (see inherited doc)
 */
ConditionItem::ItemType ConditionItem::getType() const
{
	return RunCondition;
}

void ConditionItem::updateDisplayName()
{
	std::ostringstream nameStream;
	Io::Db::ConfigVersionPtr version = runCondition->getConfigVersion();
	nameStream << version->getStructure()->getConfig()->name;
	nameStream << '/';
	nameStream << version->creationDate;
	if (not version->name.empty())
	{
		nameStream << " (" << version->name << ')';
	}
	name = nameStream.str();
}

}
}
