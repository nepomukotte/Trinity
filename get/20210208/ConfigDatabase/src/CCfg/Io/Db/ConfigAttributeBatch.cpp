#include "ConfigAttributeBatch.h"
#include <stdexcept>

#define REPORT_EXCEP throw std::logic_error(std::string("ConfigAttributeBatch: invalid data request"))

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigAttributeBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConfigAttributeBatch::ConfigAttributeBatch()
{
}

/**
 * Destructor.
 */
ConfigAttributeBatch::~ConfigAttributeBatch()
{
}

/**
 * Creates a ConfigAttribute object for a data row.
 * @param index the index of the data row.
 */
ConfigAttributePtr ConfigAttributeBatch::createEntity(size_t index)
{
	ConfigAttributePtr entity(new ConfigAttribute);
	assignEntity(index, *entity);
	return entity;
}

/**
 * Copies a data row into a ConfigAttribute object.
 * @param index the index of the data row.
 * @param entity the object where data will be copied.
 */
void ConfigAttributeBatch::assignEntity(size_t index, ConfigAttribute& entity)
{
	assert(ids.size() > index);
	assert(size() > index);
	entity.id = ids[index];
	entity.fkParent = fkParents[index];
	entity.fkVersion = fkVersions[index];
	entity.fkClassInfo = fkClassInfos[index];
	entity.index = indexes[index];
	entity.value = values[index];
}

/*
 * (see inherited doc)
 */
void ConfigAttributeBatch::append(const AbstractEntity& entity)
{
	const ConfigAttribute& cfgAttr = static_cast<const ConfigAttribute&>(entity);
	ids.push_back(cfgAttr.id);
	fkParents.push_back(cfgAttr.fkParent);
	fkVersions.push_back(cfgAttr.fkVersion);
	fkClassInfos.push_back(cfgAttr.fkClassInfo);
	indexes.push_back(cfgAttr.index);
	values.push_back(cfgAttr.value);
}

/*
 * (see inherited doc)
 */
size_t ConfigAttributeBatch::size() const
{
	assert(ids.size() == fkParents.size());
	assert(fkParents.size() == fkVersions.size());
	assert(fkVersions.size() == fkClassInfos.size());
	assert(fkClassInfos.size() == indexes.size());
	assert(indexes.size() == values.size());
	return values.size();
}

/*
 * (see inherited doc)
 */
void ConfigAttributeBatch::resize(size_t size)
{
	ids.resize(size);
	fkParents.resize(size);
	fkVersions.resize(size);
	fkClassInfos.resize(size);
	indexes.resize(size);
	values.resize(size);
}

/*
 * (see inherited doc)
 */
UidBatch& ConfigAttributeBatch::getUidBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigAttribute::ID:
		return ids;
	case ConfigAttribute::FK_PARENT:
		return fkParents;
	case ConfigAttribute::FK_VERSION:
		return fkVersions;
	case ConfigAttribute::FK_CLASSINFO:
		return fkClassInfos;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
IntBatch& ConfigAttributeBatch::getIntBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
StringBatch& ConfigAttributeBatch::getStringBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigAttribute::INDEX:
		return indexes;
	case ConfigAttribute::VALUE:
		return values;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
DateBatch& ConfigAttributeBatch::getDateBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

}
}
}
