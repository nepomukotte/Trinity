#include "SiteBatch.h"
#include <stdexcept>

#define REPORT_EXCEP throw std::logic_error(std::string("SiteBatch: invalid data request"))

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::SiteBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
SiteBatch::SiteBatch()
{
}

/**
 * Destructor.
 */
SiteBatch::~SiteBatch()
{
}

/**
 * Creates a Site object for a data row.
 * @param index the index of the data row.
 */
SitePtr SiteBatch::createEntity(size_t index)
{
	SitePtr entity(new Site);
	assignEntity(index, *entity);
	return entity;
}

/**
 * Copies a data row into a Site object.
 * @param index the index of the data row.
 * @param entity the object where data will be copied.
 */
void SiteBatch::assignEntity(size_t index, Site& entity)
{
	assert(ids.size() > index);
	assert(size() > index);
	entity.id = ids[index];
	entity.name = names[index];
}

/*
 * (see inherited doc)
 */
void SiteBatch::append(const AbstractEntity& entity)
{
	const Site& site = static_cast<const Site&>(entity);
	ids.push_back(site.id);
	names.push_back(site.name);
}

/*
 * (see inherited doc)
 */
size_t SiteBatch::size() const
{
	assert(ids.size() == names.size());
	return names.size();
}

/*
 * (see inherited doc)
 */
void SiteBatch::resize(size_t size)
{
	ids.resize(size);
	names.resize(size);
}

/*
 * (see inherited doc)
 */
UidBatch& SiteBatch::getUidBatch(int columnIdx)
{
	if(columnIdx != Site::ID)
	{
		REPORT_EXCEP;
	}
	return ids;
}

/*
 * (see inherited doc)
 */
IntBatch& SiteBatch::getIntBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
StringBatch& SiteBatch::getStringBatch(int columnIdx)
{
	if(columnIdx != Site::NAME)
	{
		REPORT_EXCEP;
	}
	return names;
}

/*
 * (see inherited doc)
 */
DateBatch& SiteBatch::getDateBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}


}
}
}
