#ifndef CCFG_IO_DB_CONFIGCLASSINFOBATCH_H_INCLUDED
#define CCFG_IO_DB_CONFIGCLASSINFOBATCH_H_INCLUDED

#include "ConfigClassInfo.h"
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
 * This class allows to efficiently transfer ConfigClassInfo objects to and from
 * the database.
 *
 * @author fchateau
 */
class ConfigClassInfoBatch : public AbstractEntityBatch
{
	UidBatch ids;
	UidBatch fkClasses;
	UidBatch fkStructures;
	IntBatch types;
	StringBatch units;
	StringBatch ranges;

public:
	ConfigClassInfoBatch();
	virtual ~ConfigClassInfoBatch();

	ConfigClassInfoPtr createEntity(size_t index);
	void assignEntity(size_t index, ConfigClassInfo& entity);

	void append(const AbstractEntity& entity);
	size_t size() const;
	void resize(size_t size);
	UidBatch& getUidBatch(int columnIdx);
	IntBatch& getIntBatch(int columnIdx);
	StringBatch& getStringBatch(int columnIdx);
	DateBatch& getDateBatch(int columnIdx);

private:
	ConfigClassInfoBatch(const ConfigClassInfoBatch& rvalue);
	ConfigClassInfoBatch& operator=(const ConfigClassInfoBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_CONFIGCLASSINFOBATCH_H_INCLUDED
