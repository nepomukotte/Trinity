#ifndef CCFG_IO_DB_CONFIGCLASSBATCH_H_INCLUDED
#define CCFG_IO_DB_CONFIGCLASSBATCH_H_INCLUDED

#include "ConfigClass.h"
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
 * This class allows to efficiently transfer ConfigClass objects to and from
 * the database.
 *
 * @author fchateau
 */
class ConfigClassBatch : public AbstractEntityBatch
{
	UidBatch ids;
	UidBatch fkParents;
	UidBatch fkCreationStructures;
	StringBatch names;
	StringBatch paths;

public:
	ConfigClassBatch();
	virtual ~ConfigClassBatch();

	ConfigClassPtr createEntity(size_t index);
	void assignEntity(size_t index, ConfigClass& entity);

	void append(const AbstractEntity& entity);
	size_t size() const;
	void resize(size_t size);
	UidBatch& getUidBatch(int columnIdx);
	IntBatch& getIntBatch(int columnIdx);
	StringBatch& getStringBatch(int columnIdx);
	DateBatch& getDateBatch(int columnIdx);

private:
	ConfigClassBatch(const ConfigClassBatch& rvalue);
	ConfigClassBatch& operator=(const ConfigClassBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_CONFIGCLASSBATCH_H_INCLUDED
