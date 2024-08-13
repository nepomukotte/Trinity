#ifndef CCFG_IO_DB_RUNBATCH_H_INCLUDED
#define CCFG_IO_DB_RUNBATCH_H_INCLUDED

#include "Run.h"
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
 * This class allows to efficiently transfer Run objects to and from the database.
 *
 * @author sizun
 */
class RunBatch : public AbstractEntityBatch
{
	UidBatch ids;
	UidBatch sites;
	UidBatch numbers;

public:
	RunBatch();
	virtual ~RunBatch();

	RunPtr createEntity(size_t index);
	void assignEntity(size_t index, Run& entity);

	void append(const AbstractEntity& entity);
	size_t size() const;
	void resize(size_t size);
	UidBatch& getUidBatch(int columnIdx);
	IntBatch& getIntBatch(int columnIdx);
	StringBatch& getStringBatch(int columnIdx);
	DateBatch& getDateBatch(int columnIdx);

private:
	RunBatch(const RunBatch& rvalue);
	RunBatch& operator=(const RunBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_RUNBATCH_H_INCLUDED
