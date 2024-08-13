#include "ConfigStructureVersionBatch.h"
#include <stdexcept>

#define REPORT_EXCEP throw std::logic_error(std::string("ConfigStructureVersionBatch: invalid data request"))

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigStructureVersionBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConfigStructureVersionBatch::ConfigStructureVersionBatch()
{
}

/**
 * Destructor.
 */
ConfigStructureVersionBatch::~ConfigStructureVersionBatch()
{
}

/**
 * Creates a ConfigStructureVersion object for a data row.
 * @param index the index of the data row.
 */
ConfigStructureVersionPtr ConfigStructureVersionBatch::createEntity(size_t index)
{
	ConfigStructureVersionPtr entity(new ConfigStructureVersion);
	assignEntity(index, *entity);
	return entity;
}

/**
 * Copies a data row into a ConfigStructureVersion object.
 * @param index the index of the data row.
 * @param entity the object where data will be copied.
 */
void ConfigStructureVersionBatch::assignEntity(size_t index, ConfigStructureVersion& entity)
{
	assert(ids.size() > index);
	assert(size() > index);
	entity.id = ids[index];
	entity.creationDate = creationDates[index];
	entity.fkConfig = fkConfigs[index];
}

/*
 * (see inherited doc)
 */
void ConfigStructureVersionBatch::append(const AbstractEntity& entity)
{
	const ConfigStructureVersion& structVer = static_cast<const ConfigStructureVersion&>(entity);
	ids.push_back(structVer.id);
	creationDates.push_back(structVer.creationDate);
	fkConfigs.push_back(structVer.fkConfig);
}

/*
 * (see inherited doc)
 */
size_t ConfigStructureVersionBatch::size() const
{
	assert(ids.size() == creationDates.size());
	assert(creationDates.size() == fkConfigs.size());
	return fkConfigs.size();
}

/*
 * (see inherited doc)
 */
void ConfigStructureVersionBatch::resize(size_t size)
{
	ids.resize(size);
	creationDates.resize(size);
	fkConfigs.resize(size);
}

/*
 * (see inherited doc)
 */
UidBatch& ConfigStructureVersionBatch::getUidBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigStructureVersion::ID:
		return ids;
	case ConfigStructureVersion::FK_CONFIG:
		return fkConfigs;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
IntBatch& ConfigStructureVersionBatch::getIntBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
StringBatch& ConfigStructureVersionBatch::getStringBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
DateBatch& ConfigStructureVersionBatch::getDateBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigStructureVersion::CREATION_DATE:
		return creationDates;
	default:
		REPORT_EXCEP;
	}
}

}
}
}
