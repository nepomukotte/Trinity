#include "RunBatch.h"
#include <stdexcept>

#define REPORT_EXCEP throw std::logic_error(std::string("RunBatch: invalid data request"))

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::RunBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
RunBatch::RunBatch()
{
}

/**
 * Destructor.
 */
RunBatch::~RunBatch()
{
}

/**
 * Creates a Run object for a data row.
 * @param index the index of the data row.
 */
RunPtr RunBatch::createEntity(size_t index)
{
	RunPtr entity(new Run);
	assignEntity(index, *entity);
	return entity;
}

/**
 * Copies a data row into a Run object.
 * @param index the index of the data row.
 * @param entity the object where data will be copied.
 */
void RunBatch::assignEntity(size_t index, Run& entity)
{
	assert(ids.size() > index);
	assert(size() > index);
	entity.id = ids[index];
	entity.fkSite = sites[index];
	entity.number = numbers[index];
}

/*
 * (see inherited doc)
 */
void RunBatch::append(const AbstractEntity& entity)
{
	const Run& run = static_cast<const Run&>(entity);
	ids.push_back(run.id);
	sites.push_back(run.fkSite);
	numbers.push_back(run.number);
}

/*
 * (see inherited doc)
 */
size_t RunBatch::size() const
{
	assert(ids.size() == numbers.size() and ids.size() == sites.size());
	return ids.size();
}

/*
 * (see inherited doc)
 */
void RunBatch::resize(size_t size)
{
	ids.resize(size);
	sites.resize(size);
	numbers.resize(size);
}

/*
 * (see inherited doc)
 */
UidBatch& RunBatch::getUidBatch(int columnIdx)
{
	if(columnIdx == Run::ID)
	{
		return ids;
	}
	else if (columnIdx == Run::FK_SITE)
	{
		return sites;
	}
	else if (columnIdx == Run::NUMBER)
	{
		return numbers;
	}

	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
IntBatch& RunBatch::getIntBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
StringBatch& RunBatch::getStringBatch(int /* columnIdx */)
{
		REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
DateBatch& RunBatch::getDateBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}


}
}
}
