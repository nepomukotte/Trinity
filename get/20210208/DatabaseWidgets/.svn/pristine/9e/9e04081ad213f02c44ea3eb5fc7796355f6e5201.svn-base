#include "VersionItem.h"
#include "ConfigItem.h"
#include "SchemaItem.h"
#include "ConnectionItem.h"
#include "Utils/Date/StdTmHelpers.h"
#include <sstream>
#include <stdexcept>

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::VersionItem class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
VersionItem::VersionItem(ConfigItem* parent, const Io::Db::ConfigVersionPtr& version, bool structColor)
	: parent(parent), version(version), name(), structureColor(structColor)
{
	updateDisplayName();
}

/**
 * Destructor.
 */
VersionItem::~VersionItem()
{
}

/*
 * (see inherited doc)
 */
bool VersionItem::isLeaf() const
{
	return true;
}

/*
 * (see inherited doc)
 */
bool VersionItem::canSetName() const
{
	return true;
}

/*
 * (see inherited doc)
 */
void VersionItem::setName(const std::string& newName)
{
	getConfig()->getSchema()->getDao().rename(*version, newName);
	updateDisplayName();
}

/*
 * (see inherited doc)
 */
const std::string& VersionItem::getName() const
{
	return version->name;
}

/*
 * (see inherited doc)
 */
const std::string& VersionItem::getDisplayName() const
{
	return name;
}

/*
 * (see inherited doc)
 */
size_t VersionItem::getItemCount()
{
	return 0;
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* VersionItem::getItem(size_t /*index*/)
{
	return 0;
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* VersionItem::getParent()
{
	return parent;
}

/*
 * (see inherited doc)
 */
void VersionItem::refresh()
{

}

/*
 * (see inherited doc)
 */
void VersionItem::removeItem(size_t /*index*/, RemoveOptions /*options*/)
{
}

/*
 * (see inherited doc)
 */
void VersionItem::clear()
{
}

/*
 * (see inherited doc)
 */
std::string VersionItem::buildUrl() const
{
	ConnectionItem* connec = parent->getSchema()->getConnection();
	return Io::DatabaseStorage::formatUrl(connec->getDbmsType(),
			connec->getUser(),
			std::string(),
			connec->getHost(),
			connec->getPort(),
			connec->getSID(),
			parent->getSchema()->getName(),
			parent->getName(),
			&version->creationDate);
}

/*
 * (see inherited doc)
 */
VersionItem::ItemType VersionItem::getType() const
{
	return Version;
}

void VersionItem::updateDisplayName()
{
	std::ostringstream nameStream;
	nameStream << version->creationDate;
	if(not version->name.empty())
	{
		nameStream << " (" << version->name << ')';
	}
	name = nameStream.str();
}

}
}
