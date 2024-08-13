#include "InteractiveLoginCallback.h"
#include "LoginDialog.h"

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::InteractiveLoginCallback class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
InteractiveLoginCallback::InteractiveLoginCallback(QWidget* parent)
	: Io::AbstractLoginCallback(), dialog(new LoginDialog(parent))
{
}

/*
 * (see inherited doc)
 */
bool InteractiveLoginCallback::requestLogin(Io::DbmsType dbms,
		const std::string& host, const std::string& port,
		const std::string& sid, std::string& user, std::string& pwd)
{
	if(not user.empty())
	{
		dialog->setUser(user);
		pwd = lookupPassword(dbms, host, port, sid, user);
		if(not pwd.empty())
		{
			return true;
		}
	}

	if(dialog->exec() == QDialog::Accepted)
	{
		user = dialog->getUser();
		pwd = dialog->getPassword();
		registerPassword(dbms, host, port, sid, user, pwd);
		dialog->clear();
		return true;
	}
	else
	{
		return false;
	}
}

/*
 * (see inherited doc)
 */
void InteractiveLoginCallback::onLoginFailed(Io::DbmsType dbms,
			const std::string& host,
			const std::string& port,
			const std::string& sid,
			const std::string& user)
{
	clearPassword(dbms, host, port, sid, user);
}

/**
 * Looks for the password associated with the specified connection infos in the cache.
 */
std::string InteractiveLoginCallback::lookupPassword(Io::DbmsType dbms,
		const std::string& host,
		const std::string& port,
		const std::string& sid,
		const std::string& user)
{
	std::string url = Io::DatabaseStorage::formatUrl(dbms, user, std::string(), host, port, sid);
	MapUrl2Passwd::const_iterator login = passwdCache.find(url);
	if(login != passwdCache.end())
	{
		return login->second;
	}
	else
	{
		return std::string();
	}
}

/**
 * Records a password associated with the specified connection infos.
 */
void InteractiveLoginCallback::registerPassword(Io::DbmsType dbms,
		const std::string& host,
		const std::string& port,
		const std::string& sid,
		const std::string& user,
		const std::string& pwd)
{
	std::string url = Io::DatabaseStorage::formatUrl(dbms, user, std::string(), host, port, sid);
	passwdCache[url] = pwd;
}

/**
 * Clears the password associated with the specified connection infos.
 */
void InteractiveLoginCallback::clearPassword(Io::DbmsType dbms,
		const std::string& host,
		const std::string& port,
		const std::string& sid,
		const std::string& user)
{
	std::string url = Io::DatabaseStorage::formatUrl(dbms, user, std::string(), host, port, sid);
	passwdCache.erase(url);
}

}
}
