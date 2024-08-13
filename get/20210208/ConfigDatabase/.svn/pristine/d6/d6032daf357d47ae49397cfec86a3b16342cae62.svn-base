#ifndef CCFG_IO_DB_CONFIGATTRIBUTEBATCH_H_INCLUDED
#define CCFG_IO_DB_CONFIGATTRIBUTEBATCH_H_INCLUDED

#include "ConfigAttribute.h"
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
 * This class allows to efficiently transfer ConfigAttribute objects to and from
 * the database.
 *
 * @author fchateau
 */
class ConfigAttributeBatch
{
	UidBatch ids;
	UidBatch fkParents;
	UidBatch fkVersions;
	UidBatch fkClassInfos;
	StringBatch indexes;
	StringBatch values;

public:
	ConfigAttributeBatch();
	virtual ~ConfigAttributeBatch();

	ConfigAttributePtr createEntity(size_t index);
	void assignEntity(size_t index, ConfigAttribute& entity);

	void append(const AbstractEntity& entity);
	size_t size() const;
	void resize(size_t size);
	UidBatch& getUidBatch(int columnIdx);
	IntBatch& getIntBatch(int columnIdx);
	StringBatch& getStringBatch(int columnIdx);
	DateBatch& getDateBatch(int columnIdx);

private:
	ConfigAttributeBatch(const ConfigAttributeBatch& rvalue);
	ConfigAttributeBatch& operator=(const ConfigAttributeBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_CONFIGATTRIBUTEBATCH_H_INCLUDED
