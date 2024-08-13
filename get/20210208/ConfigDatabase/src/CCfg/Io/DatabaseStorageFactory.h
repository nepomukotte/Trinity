#ifndef CCFG_IO_DATABASESTORAGEFACTORY_H_INCLUDED
#define CCFG_IO_DATABASESTORAGEFACTORY_H_INCLUDED

#include "CCfg/Io/AbstractStorageFactory.h"
#include "AbstractLoginCallback.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class creates database storages from URLs.
 *
 * @author fchateau
 */
class DatabaseStorageFactory : public AbstractStorageFactory
{
	LoginCallbackPtr loginCallback;
public:
	DatabaseStorageFactory();
	~DatabaseStorageFactory();

	std::auto_ptr<Storage> createStorage(const Url& url) const;

	const LoginCallbackPtr& getLoginCallback() const { return loginCallback; }
	void setLoginCallback(const LoginCallbackPtr& callback) { loginCallback = callback; }
};

}
}

#endif // CCFG_IO_DATABASESTORAGEFACTORY_H_INCLUDED
