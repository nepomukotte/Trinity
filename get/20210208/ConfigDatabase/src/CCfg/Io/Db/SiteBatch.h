#ifndef CCFG_IO_DB_SITEBATCH_H_INCLUDED
#define CCFG_IO_DB_SITEBATCH_H_INCLUDED

#include "Site.h"
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
 * This class allows to efficiently transfer Site objects to and from the database.
 *
 * @author fchateau
 */
class SiteBatch : public AbstractEntityBatch
{
	UidBatch ids;
	StringBatch names;

public:
	SiteBatch();
	virtual ~SiteBatch();

	SitePtr createEntity(size_t index);
	void assignEntity(size_t index, Site& entity);

	void append(const AbstractEntity& entity);
	size_t size() const;
	void resize(size_t size);
	UidBatch& getUidBatch(int columnIdx);
	IntBatch& getIntBatch(int columnIdx);
	StringBatch& getStringBatch(int columnIdx);
	DateBatch& getDateBatch(int columnIdx);

private:
	SiteBatch(const SiteBatch& rvalue);
	SiteBatch& operator=(const SiteBatch& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_SITEBATCH_H_INCLUDED
