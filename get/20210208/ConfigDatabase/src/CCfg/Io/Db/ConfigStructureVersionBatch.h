#ifndef CCFG_IO_DB_CONFIGSTRUCTUREVERSIONBATCH_H_INCLUDED
#define CCFG_IO_DB_CONFIGSTRUCTUREVERSIONBATCH_H_INCLUDED

#include "ConfigStructureVersion.h"
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
 * This class allows to efficiently transfer ConfigStructureVersion objects
 * to and from the database.
 *
 * @author fchateau
 */
class ConfigStructureVersionBatch : public AbstractEntityBatch
{
	UidBatch ids;
	DateBatch creationDates;
	UidBatch fkConfigs;

public:
	ConfigStructureVersionBatch();
	virtual ~ConfigStructureVersionBatch();

	ConfigStructureVersionPtr createEntity(size_t index);
	void assignEntity(size_t index, ConfigStructureVersion& entity);

	void append(const AbstractEntity& entity);
	size_t size() const;
	void resize(size_t size);
	UidBatch& getUidBatch(int columnIdx);
	IntBatch& getIntBatch(int columnIdx);
	StringBatch& getStringBatch(int columnIdx);
	DateBatch& getDateBatch(int columnIdx);

private:
	ConfigStructureVersionBatch(const ConfigStructureVersionBatch& rvalue);
	ConfigStructureVersionBatch& operator=(const ConfigStructureVersionBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_CONFIGSTRUCTUREVERSIONBATCH_H_INCLUDED
