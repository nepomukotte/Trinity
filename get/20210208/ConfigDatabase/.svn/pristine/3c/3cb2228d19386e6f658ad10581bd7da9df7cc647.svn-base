/*
 * @file RunConditionBatch.cpp
 * @author sizun
 * $Id$
 */
#include "RunConditionBatch.h"
#include <stdexcept>

#define REPORT_EXCEP throw std::logic_error(std::string("RunConditionBatch: invalid data request"))

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::RunConditionBatch class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
RunConditionBatch::RunConditionBatch()
{
}

/**
 * Destructor.
 */
RunConditionBatch::~RunConditionBatch()
{
}

/**
 * Creates a RunCondition object for a data row.
 * @param index the index of the data row.
 */
RunConditionPtr RunConditionBatch::createEntity(size_t index)
{
	RunConditionPtr entity(new RunCondition);
	assignEntity(index, *entity);
	return entity;
}

/**
 * Copies a data row into a RunCondition object.
 * @param index the index of the data row.
 * @param entity the object where data will be copied.
 */
void RunConditionBatch::assignEntity(size_t index, RunCondition& entity)
{
	assert(ids.size() > index);
	assert(size() > index);
	entity.id = ids[index];
	entity.fkRun = runs[index];
	entity.fkConfigVersion = conditions[index];
}

/*
 * (see inherited doc)
 */
void RunConditionBatch::append(const AbstractEntity& entity)
{
	const RunCondition& run = static_cast<const RunCondition&>(entity);
	ids.push_back(run.id);
	runs.push_back(run.fkRun);
	conditions.push_back(run.fkConfigVersion);
}

/*
 * (see inherited doc)
 */
size_t RunConditionBatch::size() const
{
	assert(ids.size() == runs.size() and ids.size() == conditions.size());
	return ids.size();
}

/*
 * (see inherited doc)
 */
void RunConditionBatch::resize(size_t size)
{
	ids.resize(size);
	runs.resize(size);
	conditions.resize(size);
}

/*
 * (see inherited doc)
 */
UidBatch& RunConditionBatch::getUidBatch(int columnIdx)
{
	if(columnIdx == RunCondition::ID)
	{
		return ids;
	}
	else if (columnIdx == RunCondition::FK_RUN)
	{
		return runs;
	}
	else if (columnIdx == RunCondition::FK_CONFIGVERSION)
	{
		return conditions;
	}

	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
IntBatch& RunConditionBatch::getIntBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
StringBatch& RunConditionBatch::getStringBatch(int /* columnIdx */)
{
		REPORT_EXCEP;
}

/*
 * (see inherited doc)
 */
DateBatch& RunConditionBatch::getDateBatch(int /* columnIdx */)
{
	REPORT_EXCEP;
}


}
}
}
