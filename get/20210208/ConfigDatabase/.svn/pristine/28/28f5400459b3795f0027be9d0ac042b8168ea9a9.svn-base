#ifndef CCFG_IO_DB_ABSTRACTENTITYBATCH_H_INCLUDED
#define CCFG_IO_DB_ABSTRACTENTITYBATCH_H_INCLUDED

#include "AbstractEntity.h"
#include <vector>
#include <string>
#include <ctime>

namespace CCfg
{
namespace Io
{
namespace Db
{

typedef std::vector<UID> UidBatch;
typedef std::vector<int> IntBatch;
typedef std::vector<std::string> StringBatch;
typedef std::vector<std::tm> DateBatch;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This interface describes a generic entity batch: it is necessary to make bulk
 * operations on the database.
 * It must be implemented for each entity.
 * Making bulk operation on a database requires to stream a separate array of
 * data for each column.
 * So the purpose of the implementations should be to slice objects using a
 * vector for each column in the append method, and to provide accessors to
 * these vectors.
 *
 * @author fchateau
 */
class AbstractEntityBatch
{
public:
	AbstractEntityBatch();
	virtual ~AbstractEntityBatch();

	virtual void append(const AbstractEntity& entity) = 0;
	virtual size_t size() const = 0;
	virtual void resize(size_t size) = 0;
	virtual UidBatch& getUidBatch(int columnIdx) = 0;
	virtual IntBatch& getIntBatch(int columnIdx) = 0;
	virtual StringBatch& getStringBatch(int columnIdx) = 0;
	virtual DateBatch& getDateBatch(int columnIdx) = 0;

private:
	AbstractEntityBatch(const AbstractEntityBatch& rvalue);
	AbstractEntityBatch& operator=(const AbstractEntityBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_ABSTRACTENTITYBATCH_H_INCLUDED
