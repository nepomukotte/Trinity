#include "ConfigClassBatch.h"
#include <stdexcept>

#define REPORT_EXCEP throw std::logic_error(std::string("ConfigClassBatch: invalid data request"))

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigClassBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConfigClassBatch::ConfigClassBatch()
{
}

/**
 * Destructor.
 */
ConfigClassBatch::~ConfigClassBatch()
{
}

/**
 * Creates a ConfigClass object for a data row.
 * @param index the index of the data row.
 */
ConfigClassPtr ConfigClassBatch::createEntity(size_t index)
{
	ConfigClassPtr entity(new ConfigClass);
	assignEntity(index, *entity);
	return entity;
}

/**
 * Copies a data row into a ConfigClass object.
 * @param index the index of the data row.
 * @param entity the object where data will be copied.
 */
void ConfigClassBatch::assignEntity(size_t index, ConfigClass& entity)
{
	assert(ids.size() > index);
	assert(size() > index);
	entity.id = ids[index];
	entity.fkParent = fkParents[index];
	entity.fkCreationStructure = fkCreationStructures[index];
	entity.name = names[index];
	entity.path = paths[index];
}

/*
 * (see inherited doc)
 */
void ConfigClassBatch::append(const AbstractEntity& entity)
{
	const ConfigClass& cfgClass = static_cast<const ConfigClass&>(entity);
	ids.push_back(cfgClass.id);
	fkParents.push_back(cfgClass.fkParent);
	fkCreationStructures.push_back(cfgClass.fkCreationStructure);
	names.push_back(cfgClass.name);
	paths.push_back(cfgClass.path);
}

/*
 * (see inherited doc)
 */
size_t ConfigClassBatch::size() const
{
	assert(ids.size() == fkParents.size());
	assert(fkParents.size() == fkCreationStructures.size());
	assert(fkCreationStructures.size() == names.size());
	assert(names.size() == paths.size());
	return paths.size();
}

/*
 * (see inherited doc)
 */
void ConfigClassBatch::resize(size_t size)
{
	ids.resize(size);
	fkParents.resize(size);
	fkCreationStructures.resize(size);
	names.resize(size);
	paths.resize(size);
}

/*
 * (see inherited doc)
 */
UidBatch& ConfigClassBatch::getUidBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigClass::ID:
		return ids;
	case ConfigClass::FK_PARENT:
		return fkParents;
	case ConfigClass::FK_CREATION_STRUCT:
		return fkCreationStructures;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
IntBatch& ConfigClassBatch::getIntBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
StringBatch& ConfigClassBatch::getStringBatch(int columnIdx)
{
	switch(columnIdx)
	{
	case ConfigClass::NAME:
		return names;
	case ConfigClass::PATH:
		return paths;
	default:
		REPORT_EXCEP;
	}
}

/*
 * (see inherited doc)
 */
DateBatch& ConfigClassBatch::getDateBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

}
}
}
