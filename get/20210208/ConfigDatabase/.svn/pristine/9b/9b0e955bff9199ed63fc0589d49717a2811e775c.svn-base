#ifndef CCFG_IO_DB_CONFIGENTITY_H_INCLUDED
#define CCFG_IO_DB_CONFIGENTITY_H_INCLUDED

#include "AbstractEntity.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
// forward declarations
class ConfigDao;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is the base class for entities of the object-relational mapping for configurations.
 *
 * @author fchateau
 */
class ConfigEntity : public AbstractEntity
{
protected:
	int refCount;
	ConfigDao* dao;

public:
	UID id;

	ConfigEntity();
	virtual ~ConfigEntity();

	void setDAO(ConfigDao* newDao) { dao = newDao; }
	ConfigDao* getDAO() { return dao; }
	const ConfigDao* getDAO() const { return dao; }

	friend void intrusive_ptr_add_ref(ConfigEntity* r);
	friend void intrusive_ptr_release(ConfigEntity* r);

private:
	ConfigEntity(const ConfigEntity& rvalue);
	ConfigEntity& operator=(const ConfigEntity& rvalue);
};



inline void intrusive_ptr_add_ref(ConfigEntity* r)
{
	++r->refCount;
}

inline void intrusive_ptr_release(ConfigEntity* r)
{
	if(--r->refCount == 0)
	{
		delete r;
	}
}

}
}
}

#endif // CCFG_IO_DB_CONFIGENTITY_H_INCLUDED
