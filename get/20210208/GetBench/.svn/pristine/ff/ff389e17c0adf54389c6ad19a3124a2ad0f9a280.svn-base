/**
 * @file ConfigDatabaseManager.h
 * @date 11 oct. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetBench software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#ifndef get_rc_ConfigDatabaseManager_h_INCLUDED
#define get_rc_ConfigDatabaseManager_h_INCLUDED

#include "get/rc/ConfigManager.h"
#include "CCfg/Io/AbstractLoginCallback.h"

// Forward declarations
namespace CCfg
{
namespace Io
{
	class Storage;
namespace Db
{
	class ConfigDao;
}
}
}

namespace get
{
namespace rc
{

class ConfigDatabaseManager: public ConfigManager
{
public:
	ConfigDatabaseManager(const std::string & repoUrl);
	virtual ~ConfigDatabaseManager() {}
	std::auto_ptr< CCfg::Object > getConfig(const CCfg::Io::Url & cfgUrl) const;
	class ConnectionError;
protected:
	struct InvalidLoginCallback;
	virtual std::auto_ptr< CCfg::Object > getSubConfig(const ConfigId::SubType subType, const std::string & subConfigId) const;
	virtual std::string addSubConfig(const ConfigId::SubType type, std::auto_ptr< CCfg::Object > cfg, const std::string & configId = "") const;
	virtual void getSubConfigIds(const ConfigId::SubType type, SubConfigIdSet & configIds) const;
	void initFactories();
	CCfg::Io::Db::ConfigDao & getConfigDao() const;
	std::string getConfigName(const ConfigId::SubType type) const;
	CCfg::Io::Url formatUrl(const ConfigId::SubType type, const std::string & configId) const;
	static CCfg::Io::Url obfuscateUrl(const CCfg::Io::Url & url);
	static bool dbFactoriesReady;
	CCfg::Io::Url dbUrl;	///< URL of the configuration database (without user and password).
	std::string password;	///< Database user.
	std::string user;		///< Database password.
};
//_________________________________________________________________________________________________
/** An exception thrown when trying to access to a configuration with an unknown ID. */
class ConfigDatabaseManager::ConnectionError : public ::utl::Exception
{
public:
	ConnectionError(const std::string & dbUrl, const std::string & details = "") : dbUrl(dbUrl), details(details) {}
	virtual ~ConnectionError() throw() {}
	void writeMessage(std::ostream& stream) const
	{
		stream << "Could not connect to configuration database '" << dbUrl << "'.";
		if (not details.empty())
			stream << ' ' << details;
	}
	const std::string dbUrl;
	const std::string details;
};
//_________________________________________________________________________________________________
struct ConfigDatabaseManager::InvalidLoginCallback : CCfg::Io::AbstractLoginCallback
{
	bool requestLogin(CCfg::Io::DbmsType dbms,
			const std::string& host, const std::string& port,
			const std::string& sid, std::string& user, std::string& pwd);
};
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
#endif /* get_rc_ConfigDatabaseManager_h_INCLUDED */
