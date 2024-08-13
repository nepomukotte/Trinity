#include "AbstractLoginCallback.h"

namespace CCfg
{
namespace Io
{

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::AbstractLoginCallback class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * @fn bool AbstractLoginCallback::requestLogin(DbmsType dbms,
 * 		const std::string& host, const std::string& port,
 *  	const std::string& sid, std::string& user, std::string& pwd)
 * This method is called each time login informations are missing to establish a
 * database connection.
 * It can be implemented to prompt the user for these informations with a dialog
 * box.
 * @param dbms the database system the program tries to connect to.
 * @param host the database server host name.
 * @param port the database server port.
 * @param sid the SID (oracle) where we connect.
 * @param user inout the user account to use. It must be defined before returning.
 * It can be already set if this information was given in the URL. In this case
 * only the password is needed, nevertheless it should be still possible to
 * change the user.
 * @param pwd out the password to use. It must be defined before returning. It
 * always empty when the function is called, otherwise connection would have been
 * made without calling this method.
 * @return Returns true to confirm login and false to cancel it.
 */

/**
 * This method is called each time the login informations provided by requestLogin
 * lead to a connection failure.
 */
void AbstractLoginCallback::onLoginFailed(DbmsType /*dbms*/,
			const std::string& /*host*/,
			const std::string& /*port*/,
			const std::string& /*sid*/,
			const std::string& /*user*/)
{
}


}
}
