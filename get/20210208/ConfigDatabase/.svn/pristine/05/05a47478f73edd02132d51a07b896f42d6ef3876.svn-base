#ifndef CCFG_IO_DB_CONFIGVERSIONBATCH_H_INCLUDED
#define CCFG_IO_DB_CONFIGVERSIONBATCH_H_INCLUDED

#include "ConfigVersion.h"
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
 * This class allows to efficiently transfer ConfigVersion objects to and from
 * the database.
 *
 * @author fchateau
 */
class ConfigVersionBatch : public AbstractEntityBatch
{
	UidBatch ids;
	DateBatch creationDates;
	StringBatch names;
	UidBatch fkStructures;

public:
	ConfigVersionBatch();
	virtual ~ConfigVersionBatch();

	ConfigVersionPtr createEntity(size_t index);
	void assignEntity(size_t index, ConfigVersion& entity);

	void append(const AbstractEntity& entity);
	size_t size() const;
	void resize(size_t size);
	UidBatch& getUidBatch(int columnIdx);
	IntBatch& getIntBatch(int columnIdx);
	StringBatch& getStringBatch(int columnIdx);
	DateBatch& getDateBatch(int columnIdx);

private:
	ConfigVersionBatch(const ConfigVersionBatch& rvalue);
	ConfigVersionBatch& operator=(const ConfigVersionBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_CONFIGVERSIONBATCH_H_INCLUDED
