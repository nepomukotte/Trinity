#include "ConfigBatch.h"
#include <stdexcept>

#define REPORT_EXCEP throw std::logic_error(std::string("ConfigBatch: invalid data request"))

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ConfigBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConfigBatch::ConfigBatch()
{
}

/**
 * Destructor.
 */
ConfigBatch::~ConfigBatch()
{
}

/**
 * Creates a configuration object for a data row.
 * @param index the index of the data row.
 */
ConfigPtr ConfigBatch::createEntity(size_t index)
{
	ConfigPtr entity(new Config);
	assignEntity(index, *entity);
	return entity;
}

/**
 * Copies a data row into a Config object.
 * @param index the index of the data row.
 * @param entity the object where data will be copied.
 */
void ConfigBatch::assignEntity(size_t index, Config& entity)
{
	assert(ids.size() > index);
	assert(size() > index);
	entity.id = ids[index];
	entity.name = names[index];
}

/*
 * (see inherited doc)
 */
void ConfigBatch::append(const AbstractEntity& entity)
{
	const Config& config = static_cast<const Config&>(entity);
	ids.push_back(config.id);
	names.push_back(config.name);
}

/*
 * (see inherited doc)
 */
size_t ConfigBatch::size() const
{
	assert(ids.size() == names.size());
	return names.size();
}

/*
 * (see inherited doc)
 */
void ConfigBatch::resize(size_t size)
{
	ids.resize(size);
	names.resize(size);
}

/*
 * (see inherited doc)
 */
UidBatch& ConfigBatch::getUidBatch(int columnIdx)
{
	if(columnIdx != Config::ID)
	{
		REPORT_EXCEP;
	}
	return ids;
}

/*
 * (see inherited doc)
 */
IntBatch& ConfigBatch::getIntBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
StringBatch& ConfigBatch::getStringBatch(int columnIdx)
{
	if(columnIdx != Config::NAME)
	{
		REPORT_EXCEP;
	}
	return names;
}

/*
 * (see inherited doc)
 */
DateBatch& ConfigBatch::getDateBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}


}
}
}
