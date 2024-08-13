#ifndef CCFG_IO_DB_SQLITE3CONFIGDAO_H_INCLUDED
#define CCFG_IO_DB_SQLITE3CONFIGDAO_H_INCLUDED

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
 * This DAO subclass redefines queries and methods specific to SQLite3 databases.
 *
 * @author psizun
 */
class Sqlite3ConfigDao : public ConfigDao
{
public:
	Sqlite3ConfigDao();
	virtual ~Sqlite3ConfigDao();
	// @{ Connection
	void connect(const std::string& dbPath);
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
	// @{ Primary key generation
	UID allocIdRange(const char* tableName, size_t count);
	void allocIdRanges(const UID (&counts)[TableCount], UID (&ids)[TableCount]);
	// @}

protected:
	std::auto_ptr<TableBuilder> createTableBuilder();

private:
	Sqlite3ConfigDao(const Sqlite3ConfigDao& rvalue);
	Sqlite3ConfigDao& operator=(const Sqlite3ConfigDao& rvalue);
};

}
}
}

#endif // CCFG_IO_DB_SQLITE3CONFIGDAO_H_INCLUDED
