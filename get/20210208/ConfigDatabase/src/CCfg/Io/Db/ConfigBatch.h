#ifndef CCFG_IO_DB_CONFIGBATCH_H_INCLUDED
#define CCFG_IO_DB_CONFIGBATCH_H_INCLUDED

#include "Config.h"
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
 * This class allows to efficiently transfer Config objects to and from the database.
 *
 * @author fchateau
 */
class ConfigBatch : public AbstractEntityBatch
{
	UidBatch ids;
	StringBatch names;

public:
	ConfigBatch();
	virtual ~ConfigBatch();

	ConfigPtr createEntity(size_t index);
	void assignEntity(size_t index, Config& entity);

	void append(const AbstractEntity& entity);
	size_t size() const;
	void resize(size_t size);
	UidBatch& getUidBatch(int columnIdx);
	IntBatch& getIntBatch(int columnIdx);
	StringBatch& getStringBatch(int columnIdx);
	DateBatch& getDateBatch(int columnIdx);

private:
	ConfigBatch(const ConfigBatch& rvalue);
	ConfigBatch& operator=(const ConfigBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_CONFIGBATCH_H_INCLUDED
