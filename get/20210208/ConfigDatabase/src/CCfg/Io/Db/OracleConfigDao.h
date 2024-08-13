#ifndef CCFG_IO_DB_ORACLECONFIGDAO_H_INCLUDED
#define CCFG_IO_DB_ORACLECONFIGDAO_H_INCLUDED

#include "ConfigDao.h"

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
 * This DAO subclass redefines queries and methods specific to Oracle databases.
 *
 * @author fchateau
 */
class OracleConfigDao : public ConfigDao
{
public:
	OracleConfigDao();
	virtual ~OracleConfigDao();
	// @{ Connection
	void connect(const std::string& connectionString);
	void connect(const std::string& host, const std::string& port, const std::string& sid, const std::string& user, const std::string& pwd);
	// @}
	// @{ Database and user management
	ConfigDatabasePtr searchDatabase(const std::string& dbName);
	ConfigDatabasePtr getCurrentDatabase();
	ConfigDatabasePtrSetPtr getDatabases();
	std::vector<std::string> getUserNames();
	void createDatabaseAndUser(const std::string& name, const std::string& pwd);
	void removeDatabaseAndUser(const std::string& name);
	void changeUserPassword(const std::string& userName, const std::string& newPassword);
	// @}
	// @{ Privileges management
	int getPrivileges(const std::string& userName, const std::string& tableName);
	void revokePrivileges(const std::string& userName);
	// @}
	// @{ Roles management
	ConfigUserRole getRole(const std::string& userName, int tableIndex);
	bool isSchemaAdmin(const std::string& userName);
	// @}
	// @{ Tables management
	bool checkTableExist(const char* tableName);
	// @}

protected:
	std::auto_ptr<TableBuilder> createTableBuilder();

private:
	OracleConfigDao(const OracleConfigDao& rvalue);
	OracleConfigDao& operator=(const OracleConfigDao& rvalue);
};


}
}
}

#endif // CCFG_IO_DB_ORACLECONFIGDAO_H_INCLUDED
