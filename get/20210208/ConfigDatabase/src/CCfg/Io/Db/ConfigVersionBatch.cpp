#include "ConfigVersionBatch.h"
#include <stdexcept>

#define REPORT_EXCEP throw std::logic_error(std::string("ConfigVersionBatch: invalid data request"))

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigVersionBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConfigVersionBatch::ConfigVersionBatch()
{
}

/**
 * Destructor.
 */
ConfigVersionBatch::~ConfigVersionBatch()
{
}

/**
 * Creates a configuration object for a data row.
 * @param index the index of the data row.
 */
ConfigVersionPtr ConfigVersionBatch::createEntity(size_t index)
{
	ConfigVersionPtr entity(new ConfigVersion);
	assignEntity(index, *entity);
	return entity;
}

/**
 * Copies a data row into a ConfigVersion object.
 * @param index the index of the data row.
 * @param entity the object where data will be copied.
 */
void ConfigVersionBatch::assignEntity(size_t index, ConfigVersion& entity)
{
	assert(ids.size() > index);
	assert(size() > index);
	entity.id = ids[index];
	entity.creationDate = creationDates[index];
	entity.name = names[index];
	entity.fkStructure = fkStructures[index];
}

/*
 * (see inherited doc)
 */
void ConfigVersionBatch::append(const AbstractEntity& entity)
{
	const ConfigVersion& version = static_cast<const ConfigVersion&>(entity);
	ids.push_back(version.id);
	creationDates.push_back(version.creationDate);
	names.push_back(version.name);
	fkStructures.push_back(version.fkStructure);
}

/*
 * (see inherited doc)
 */
size_t ConfigVersionBatch::size() const
{
	assert(ids.size() == creationDates.size());
	assert(creationDates.size() == names.size());
	assert(names.size() == fkStructures.size());
	return fkStructures.size();
}

/*
 * (see inherited doc)
 */
void ConfigVersionBatch::resize(size_t size)
{
	ids.resize(size);
	creationDates.resize(size);
	names.resize(size);
	fkStructures.resize(size);
}

/*
 * (see inherited doc)
 */
UidBatch& ConfigVersionBatch::getUidBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigVersion::ID:
		return ids;
	case ConfigVersion::FK_STRUCTURE:
		return fkStructures;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
IntBatch& ConfigVersionBatch::getIntBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
StringBatch& ConfigVersionBatch::getStringBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigVersion::NAME:
		return names;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
DateBatch& ConfigVersionBatch::getDateBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigVersion::CREATION_DATE:
		return creationDates;
	default:
		REPORT_EXCEP;
	}
}

}
}
}
