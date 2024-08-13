#ifndef CCFG_IO_ABSTRACTLOGINCALLBACK_H_INCLUDED
#define CCFG_IO_ABSTRACTLOGINCALLBACK_H_INCLUDED

#include "DatabaseStorage.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This interface is used by the database access layer to interactively prompt
 * login informations to an entity (user, program, ...).
 * This feature allows to keep the password secret by avoiding the obligation to
 * put it in URLs.
 *
 * @author fchateau
 */
class AbstractLoginCallback
{
public:
	virtual ~AbstractLoginCallback() {}
	virtual bool requestLogin(DbmsType dbms,
			const std::string& host,
			const std::string& port,
			const std::string& sid,
			std::string& user,
			std::string& pwd) = 0;

	virtual void onLoginFailed(DbmsType dbms,
			const std::string& host,
			const std::string& port,
			const std::string& sid,
			const std::string& user);
};

}
}

#endif // CCFG_IO_ABSTRACTLOGINCALLBACK_H_INCLUDED
