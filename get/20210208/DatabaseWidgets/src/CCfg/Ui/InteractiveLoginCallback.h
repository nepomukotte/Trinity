#ifndef CCFG_UI_INTERACTIVELOGINCALLBACK_H_INCLUDED
#define CCFG_UI_INTERACTIVELOGINCALLBACK_H_INCLUDED

#include "CCfg/Io/AbstractLoginCallback.h"
#include <QWidget>
#include <map>

namespace CCfg
{
namespace Ui
{
// forward declarations
class LoginDialog;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class requests login informations through a dialog box when the database
 * module request it.
 *
 * @author fchateau
 */
class InteractiveLoginCallback : public Io::AbstractLoginCallback
{
	typedef std::map<std::string, std::string> MapUrl2Passwd;
	LoginDialog* dialog;
	MapUrl2Passwd passwdCache;

public:
	InteractiveLoginCallback(QWidget* parent);
	// @{ Interface implementation
	bool requestLogin(Io::DbmsType dbms,
			const std::string& host,
			const std::string& port,
			const std::string& sid,
			std::string& user,
			std::string& pwd);
	void onLoginFailed(Io::DbmsType dbms,
				const std::string& host,
				const std::string& port,
				const std::string& sid,
				const std::string& user);
	// @}
	// @{ Password cache management
	std::string lookupPassword(Io::DbmsType dbms,
			const std::string& host,
			const std::string& port,
			const std::string& sid,
			const std::string& user);
	void registerPassword(Io::DbmsType dbms,
			const std::string& host,
			const std::string& port,
			const std::string& sid,
			const std::string& user,
			const std::string& pwd);
	void clearPassword(Io::DbmsType dbms,
			const std::string& host,
			const std::string& port,
			const std::string& sid,
			const std::string& user);
	// @}

private:
	InteractiveLoginCallback(const InteractiveLoginCallback& rvalue);
	InteractiveLoginCallback& operator=(const InteractiveLoginCallback& rvalue);
};

}
}

#endif // CCFG_UI_INTERACTIVELOGINCALLBACK_H_INCLUDED
