/*
 * @file RunConditionBatch.h
 * @author sizun
 * $Id$
 */
#ifndef CCFG_IO_DB_RUNCONDITIONBATCH_H_INCLUDED
#define CCFG_IO_DB_RUNCONDITIONBATCH_H_INCLUDED

#include "RunCondition.h"
#include "AbstractEntityBatch.h"

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
 * This class allows to efficiently transfer RunCondition objects to and from the database.
 *
 * @author sizun
 */
class RunConditionBatch : public AbstractEntityBatch
{
	UidBatch ids;
	UidBatch runs;
	UidBatch conditions;

public:
	RunConditionBatch();
	virtual ~RunConditionBatch();

	RunConditionPtr createEntity(size_t index);
	void assignEntity(size_t index, RunCondition& entity);

	void append(const AbstractEntity& entity);
	size_t size() const;
	void resize(size_t size);
	UidBatch& getUidBatch(int columnIdx);
	IntBatch& getIntBatch(int columnIdx);
	StringBatch& getStringBatch(int columnIdx);
	DateBatch& getDateBatch(int columnIdx);

private:
	RunConditionBatch(const RunConditionBatch& rvalue);
	RunConditionBatch& operator=(const RunConditionBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_RUNCONDITIONBATCH_H_INCLUDED
