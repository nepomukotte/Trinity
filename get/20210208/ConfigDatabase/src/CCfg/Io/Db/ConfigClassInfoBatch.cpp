#include "ConfigClassInfoBatch.h"
#include <stdexcept>

#define REPORT_EXCEP throw std::logic_error(std::string("ConfigClassInfoBatch: invalid data request"))

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigClassInfoBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConfigClassInfoBatch::ConfigClassInfoBatch()
{
}

/**
 * Destructor.
 */
ConfigClassInfoBatch::~ConfigClassInfoBatch()
{
}

/**
 * Creates a configuration object for a data row.
 * @param index the index of the data row.
 */
ConfigClassInfoPtr ConfigClassInfoBatch::createEntity(size_t index)
{
	ConfigClassInfoPtr entity(new ConfigClassInfo);
	assignEntity(index, *entity);
	return entity;
}

/**
 * Copies a data row into a ConfigClassInfo object.
 * @param index the index of the data row.
 * @param entity the object where data will be copied.
 */
void ConfigClassInfoBatch::assignEntity(size_t index, ConfigClassInfo& entity)
{
	assert(ids.size() > index);
	assert(size() > index);
	entity.id = ids[index];
	entity.fkClass = fkClasses[index];
	entity.fkStructure = fkStructures[index];
	entity.type = static_cast<AttributeTypeId>(types[index]);
	entity.unit = units[index];
	entity.range = ranges[index];
}

/*
 * (see inherited doc)
 */
void ConfigClassInfoBatch::append(const AbstractEntity& entity)
{
	const ConfigClassInfo& classInfo = static_cast<const ConfigClassInfo&>(entity);
	ids.push_back(classInfo.id);
	fkClasses.push_back(classInfo.fkClass);
	fkStructures.push_back(classInfo.fkStructure);
	types.push_back(classInfo.type);
	units.push_back(classInfo.unit);
	ranges.push_back(classInfo.range);
}

/*
 * (see inherited doc)
 */
size_t ConfigClassInfoBatch::size() const
{
	assert(ids.size() == fkClasses.size());
	assert(fkClasses.size() == fkStructures.size());
	assert(fkStructures.size() == types.size());
	assert(types.size() == units.size());
	assert(units.size() == ranges.size());
	return ranges.size();
}

/*
 * (see inherited doc)
 */
void ConfigClassInfoBatch::resize(size_t size)
{
	ids.resize(size);
	fkClasses.resize(size);
	fkStructures.resize(size);
	types.resize(size);
	units.resize(size);
	ranges.resize(size);
}

/*
 * (see inherited doc)
 */
UidBatch& ConfigClassInfoBatch::getUidBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigClassInfo::ID:
		return ids;
	case ConfigClassInfo::FK_CLASS:
		return fkClasses;
	case ConfigClassInfo::FK_STRUCTURE:
		return fkStructures;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
IntBatch& ConfigClassInfoBatch::getIntBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigClassInfo::TYPE:
		return types;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
StringBatch& ConfigClassInfoBatch::getStringBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigClassInfo::UNIT:
		return units;
	case ConfigClassInfo::RANGE:
		return ranges;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
DateBatch& ConfigClassInfoBatch::getDateBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

}
}
}
