#ifndef CCFG_IO_DB_ABSTRACTENTITY_H_INCLUDED
#define CCFG_IO_DB_ABSTRACTENTITY_H_INCLUDED

namespace CCfg
{
namespace Io
{
namespace Db
{
typedef unsigned long UID;
extern const UID NULL_UID;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is an interface that every entity mapping a database table should implement.
 * @author fchateau
 */
class AbstractEntity
{
public:
	AbstractEntity();
	virtual ~AbstractEntity();
	virtual UID getId() const = 0;

private:
	AbstractEntity(const AbstractEntity& rvalue);
	AbstractEntity& operator=(const AbstractEntity& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_ABSTRACTENTITY_H_INCLUDED
