#include "OracleConfigDao.h"
#include "Transaction.h"
#include "OracleTableBuilder.h"
#include "ObjRelMapping.hpp"
#include "Utils/String/Comparison.hpp"
#include <oracle/soci-oracle.h>
#include <soci.h>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <cctype>

namespace CCfg
{
namespace Io
{
namespace Db
{
using namespace SOCI;

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::OracleConfigDao class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
OracleConfigDao::OracleConfigDao()
{
}

/**
 * Destructor.
 */
OracleConfigDao::~OracleConfigDao()
{
}

/*
 * (see inherited doc)
 */
void OracleConfigDao::connect(const std::string& connectionString)
{
	session.reset(new Session(oracle, connectionString));
	session->once << "SELECT user FROM dual", into(dbName);
}

/**
 * Connects to an Oracle database (using RDMS specific parameters).
 */
void OracleConfigDao::connect(const std::string& host, const std::string& port,
		const std::string& sid, const std::string& user, const std::string& pwd)
{
	std::ostringstream connectionString;
	connectionString << "service=" << host << ':' << port << '/' << sid
		<< " user=" << user
		<< " password=" << pwd;
	if(Utils::String::EqualNoCaseStr()(user, "SYS"))
	{
		connectionString << " mode=sysdba";
	}
	connect(connectionString.str());
	dbName = user;
}

/*
 * (see inherited doc)
 */
ConfigDatabasePtr OracleConfigDao::searchDatabase(const std::string& dbName)
{
	Statement stmt = (session->prepare <<
		"SELECT username AS schema_name "
		"FROM ALL_USERS "
		"WHERE username=:name",
		use(const_cast<std::string&>(dbName)));
	return loadDatabase(stmt);
}

/*
 * (see inherited doc)
 */
ConfigDatabasePtr OracleConfigDao::getCurrentDatabase()
{
	ConfigDatabasePtr dbObj(new ConfigDatabase);
	dbObj->setDAO(this);
	session->once << "SELECT user FROM dual", into(dbObj->name);
	return dbObj;
}

/*
 * (see inherited doc)
 */
ConfigDatabasePtrSetPtr OracleConfigDao::getDatabases()
{
	ConfigDatabasePtrSetPtr databases(new ConfigDatabasePtrSet);
	Rowset<ConfigDatabasePtr> rowset =
		(session->prepare << "SELECT username AS schema_name FROM ALL_USERS");
	for(Rowset<ConfigDatabasePtr>::iterator it = rowset.begin(); it != rowset.end(); ++it)
	{
		(*it)->setDAO(this);
		databases->insert(*it);
	}
	return databases;
}

/*
 * (see inherited doc)
 */
std::vector<std::string> OracleConfigDao::getUserNames()
{
	std::vector<std::string> users;
	Rowset<std::string> rowset = (session->prepare <<
			"SELECT username FROM ALL_USERS");
	for(Rowset<std::string>::iterator it = rowset.begin(); it != rowset.end(); ++it)
	{
		users.push_back(*it);
	}
	return users;
}

/*
 * (see inherited doc)
 */
void OracleConfigDao::createDatabaseAndUser(const std::string& name, const std::string& pwd)
{
	session->once << "CREATE USER " << name << " IDENTIFIED BY " << pwd;
	session->once << "GRANT CONNECT, RESOURCE TO " << name;
}

/*
 * (see inherited doc)
 */
void OracleConfigDao::removeDatabaseAndUser(const std::string& name)
{
	session->once << "DROP USER " << name << " CASCADE";
}

/*
 * (see inherited doc)
 */
void OracleConfigDao::changeUserPassword(const std::string& userName, const std::string& newPassword)
{
	session->once <<
		"ALTER USER " << userName << " "
		"IDENTIFIED BY " << newPassword;
}

/*
 * (see inherited doc)
 */
void OracleConfigDao::revokePrivileges(const std::string& userName)
{
	for(size_t i = 0; i < TableCount+1; ++i)
	{
		session->once <<
			"REVOKE ALL "
			"ON " << dbName << '.' << TABLES[i] << " "
			"FROM " << userName;
	}
}

/*
 * (see inherited doc)
 */
int OracleConfigDao::getPrivileges(const std::string& userName, const std::string& tableName)
{
	eIndicator ind;
	int privileges = 0;
	ConfigDatabasePtr user = getCurrentDatabase();
	if(user->name == "SYS")
	{
		session->once <<
			"SELECT SUM(CASE PRIVILEGE "
				"WHEN 'SELECT' THEN 1 "
				"WHEN 'INSERT' THEN 2 "
				"WHEN 'UPDATE' THEN 4 "
				"WHEN 'DELETE' THEN 8 "
				"ELSE 0 END) AS AccessFlag "
			"FROM DBA_TAB_PRIVS "
			"WHERE table_name=UPPER(:tblName) AND owner=:db AND grantee=:usr",
			use(const_cast<std::string&>(tableName)),
			use(const_cast<std::string&>(dbName)),
			use(const_cast<std::string&>(userName)),
			into(privileges, ind);
	}
	else
	{
		session->once <<
			"SELECT SUM(CASE PRIVILEGE "
				"WHEN 'SELECT' THEN 1 "
				"WHEN 'INSERT' THEN 2 "
				"WHEN 'UPDATE' THEN 4 "
				"WHEN 'DELETE' THEN 8 "
				"ELSE 0 END) AS AccessFlag "
			"FROM ALL_TAB_PRIVS "
			"WHERE table_name=UPPER(:tblName) AND table_schema=:db AND grantee=:usr",
			use(const_cast<std::string&>(tableName)),
			use(const_cast<std::string&>(dbName)),
			use(const_cast<std::string&>(userName)),
			into(privileges, ind);
	}

	return privileges;
}

/*
 * (see inherited doc)
 */
bool OracleConfigDao::isSchemaAdmin(const std::string& userName)
{
	return userName == dbName;
}

/*
 * (see inherited doc)
 */
bool OracleConfigDao::checkTableExist(const char* tableName)
{
	int exist = -1;
	std::string name;
	size_t len = strlen(tableName);
	name.reserve(len);
	std::transform(&tableName[0], &tableName[len], std::back_inserter(name),
			static_cast<int (*)(int)>(toupper));

	session->once << "SELECT COUNT(*) AS tableExists "
		"FROM ALL_TABLES "
		"WHERE OWNER=:schema AND TABLE_NAME=:name",
		into(exist),
		use(dbName),
		use(name);
	return exist == 1;
}

/*
 * (see inherited doc)
 */
std::auto_ptr<TableBuilder> OracleConfigDao::createTableBuilder()
{
	return std::auto_ptr<TableBuilder>(new OracleTableBuilder);
}

}
}
}
