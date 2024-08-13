#include "PrimaryKeyGenerator.h"
#include <algorithm>
#include <functional>
#include <stdexcept>

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::PrimaryKeyGenerator class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
PrimaryKeyGenerator::PrimaryKeyGenerator(ConfigDao* dao) : dao(dao)
{
}

/**
 * Destructor.
 */
PrimaryKeyGenerator::~PrimaryKeyGenerator()
{
}

/**
 * Reserves a primary key range for a single table.
 * @param range out the range that will be initialized.
 * @param table the table for which the range should be reserved.
 * @param count the number of elements to reserve for the range.
 */
void PrimaryKeyGenerator::reserve(SingleRange& range, TableId table, UID count)
{
	switch(table)
	{
	case TConfig:range.current = dao->allocIdRangeForConfigs(count);
		break;
	case TStructureVersion:range.current = dao->allocIdRangeForStructureVersions(count);
		break;
	case TVersion:range.current = dao->allocIdRangeForVersions(count);
		break;
	case TClass:range.current = dao->allocIdRangeForClasses(count);
		break;
	case TClassInfo:range.current = dao->allocIdRangeForClassInfos(count);
		break;
	case TAttribute:range.current = dao->allocIdRangeForAttributes(count);
		break;
	default:
		throw std::logic_error("TableId out of range");
	}
	range.limit = range.current + count;
}

/**
 * Reserves a primary key range for all tables.
 * @param range out the range that will be initialized.
 * @param counts the number of elements to reserve in each table for the range.
 */
void PrimaryKeyGenerator::reserve(MultiRange& range, UID (&counts)[TableCount])
{
	dao->allocIdRanges(range.current, counts);
	// vectorial form of: range.limit = range.current + counts
	std::transform(range.current, range.current + TableCount,
		counts, range.limit,
	    std::plus<UID>());
}

/**
 * Reserves a primary key range for all tables.
 */
void PrimaryKeyGenerator::reserve(MultiRange& range, UID nConfig, UID nStruct, UID nVersion, UID nClass, UID nInfos, UID nAttr)
{
	UID counts[TableCount] = {nConfig, nStruct, nVersion, nClass, nInfos, nAttr};
	dao->allocIdRanges(range.current, counts);
	// vectorial form of: range.limit = range.current + counts
	std::transform(range.current, range.current + TableCount,
		counts, range.limit,
	    std::plus<UID>());
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::PrimaryKeyGenerator::SingleRange class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
PrimaryKeyGenerator::SingleRange::SingleRange()
	: current(NULL_UID), limit(NULL_UID)
{
}

/**
 * Full constructor.
 */
PrimaryKeyGenerator::SingleRange::SingleRange(UID poolBegin, UID poolSize)
	: current(poolBegin), limit(poolBegin + poolSize)
{
}

/**
 * Returns the next primary key value in the range.
 * This function should be called each time a new key is need, even if the range
 * has just been initialized (because currentVal is not valid yet).
 * @throws PrimaryKeyGenerator::RangeDepleted if the number of request values
 * is now greater than the number of reserved values.
 */
UID PrimaryKeyGenerator::SingleRange::nextVal()
{
	if(not hasNext())
	{
		throw RangeDepleted();
	}
	return ++current;
}

/**
 * Returns the primary key returned by the lastest call to nextVal.
 */
UID PrimaryKeyGenerator::SingleRange::currentVal() const
{
	return current;
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::PrimaryKeyGenerator::MultiRange class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
PrimaryKeyGenerator::MultiRange::MultiRange()
	: current(), limit()
{
	std::fill(current, current + TableCount, NULL_UID);
	std::fill(current, current + TableCount, NULL_UID);
}

/**
 * Full constructor.
 */
PrimaryKeyGenerator::MultiRange::MultiRange(UID (&poolBegin)[TableCount], UID (&poolSizes)[TableCount])
	: current(), limit()
{
	for(size_t i = 0; i < TableCount; ++i)
	{
		current[i] = poolBegin[i];
		limit[i] = poolBegin[i] + poolSizes[i];
	}
}

/**
 * Returns the next primary key value in the range for a table.
 * This function should be called each time a new key is need, even if the range
 * has just been initialized (because currentVal is not valid yet).
 * @throws PrimaryKeyGenerator::RangeDepleted if the number of request values
 * for the table is now greater than the number of reserved values.
 */
UID PrimaryKeyGenerator::MultiRange::nextVal(TableId table)
{
	if(not hasNext(table))
	{
		throw RangeDepleted();
	}
	return ++current[table];
}

/**
 * Returns the primary key returned by the lastest call to nextVal for a table.
 */
UID PrimaryKeyGenerator::MultiRange::currentVal(TableId table) const
{
	return current[table];
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db:: class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
PrimaryKeyGenerator::RangeDepleted::RangeDepleted() : Exception()
{
}

/**
 * Copy constructor.
 */
PrimaryKeyGenerator::RangeDepleted::RangeDepleted(const RangeDepleted& rvalue)
	: Exception(rvalue)
{
}

/**
 * Assignment operator;
 */
PrimaryKeyGenerator::RangeDepleted&
PrimaryKeyGenerator::RangeDepleted::operator=(const RangeDepleted& rvalue)
{
	Exception::operator=(rvalue);
	return *this;
}

/*
 * (see inherited doc)
 */
std::string
PrimaryKeyGenerator::RangeDepleted::formatMessage() const
{
	return "Cannot fetch the next number from the primary key range because it is depleted";
}

}
}
}
