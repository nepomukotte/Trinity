#ifndef CCFG_IO_DB_PRIMARYKEYGENERATOR_H_INCLUDED
#define CCFG_IO_DB_PRIMARYKEYGENERATOR_H_INCLUDED

#include "ConfigDao.h"
#include "Exception.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * Allocates ranges of primary keys from the database.
 *
 * @author fchateau
 */
class PrimaryKeyGenerator
{
	ConfigDao* dao;
public:
	// inner classes
	class SingleRange;
	class MultiRange;
	class RangeDepleted;

	PrimaryKeyGenerator(ConfigDao* dao);
	~PrimaryKeyGenerator();

	void reserve(SingleRange& range, TableId table, UID count);
	void reserve(MultiRange& range, UID (&counts)[TableCount]);
	void reserve(MultiRange& range, UID nConfig, UID nStruct, UID nVersion, UID nClass, UID nInfos, UID nAttr);
};

/**
 * Represents a range of primary keys allocated from the database for a single table.
 */
class PrimaryKeyGenerator::SingleRange
{
	friend class PrimaryKeyGenerator;
	UID current;
	UID limit;
public:
	SingleRange();
	SingleRange(UID poolBegin, UID poolSize);
	~SingleRange() {}
	UID nextVal();
	UID currentVal() const;
	bool hasNext() const { return current < limit; }

private:
	SingleRange(const SingleRange& rvalue);
	SingleRange& operator=(const SingleRange& rvalue);
};

/**
 * Represents a range of primary keys allocated from the database for all tables.
 */
class PrimaryKeyGenerator::MultiRange
{
	friend class PrimaryKeyGenerator;
	UID current[TableCount];
	UID limit[TableCount];
public:
	MultiRange();
	MultiRange(UID (&poolBegin)[TableCount], UID (&poolSizes)[TableCount]);
	~MultiRange() {}
	UID nextVal(TableId table);
	UID currentVal(TableId table) const;
	bool hasNext(TableId table) const { return current[table] < limit[table]; }

private:
	MultiRange(const MultiRange& rvalue);
	MultiRange& operator=(const MultiRange& rvalue);
};

/**
 * Reports that a range has run out of keys to return.
 */
class PrimaryKeyGenerator::RangeDepleted : public Exception
{
public:
	RangeDepleted();
	RangeDepleted(const RangeDepleted& rvalue);
	~RangeDepleted() throw() {}
	RangeDepleted& operator=(const RangeDepleted& rvalue);

private:
	std::string formatMessage() const;
};

}
}
}

#endif // CCFG_IO_DB_PRIMARYKEYGENERATOR_H_INCLUDED
