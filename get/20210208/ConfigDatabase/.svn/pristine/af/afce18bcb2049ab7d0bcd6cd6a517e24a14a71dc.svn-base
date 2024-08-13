#include "MySqlConfigDao.h"
#include "Transaction.h"
#include "MySqlTableBuilder.h"
#include "ObjRelMapping.hpp"
#include <mysql/soci-mysql.h>
#include <soci.h>
#include <sstream>

namespace CCfg
{
namespace Io
{
namespace Db
{
using namespace SOCI;

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::MySqlConfigDao class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
MySqlConfigDao::MySqlConfigDao() : ConfigDao()
{
}

/**
 * Destructor.
 */
MySqlConfigDao::~MySqlConfigDao()
{
}

/*
 * (see inherited doc)
 */
void MySqlConfigDao::connect(const std::string& connectionString)
{
	session.reset(new Session(mysql, connectionString));
	eIndicator ind;
	session->once << "SELECT database() AS currentDb", into(dbName, ind);
}

/**
 * Connects to a MySQL database (using RDMS specific parameters).
 */
void MySqlConfigDao::connect(const std::string& host, const std::string& port, const std::string& uid, const std::string& pwd, const std::string& db)
{
	std::ostringstream connectionString;
	connectionString << "host=" << host;
	if(not port.empty())
	{
		connectionString << " port=" << port;
	}
	connectionString << " user=" << uid
		<< " password=" << pwd;
	if(not db.empty())
	{
		connectionString << " db=" << db;
	}
	connect(connectionString.str());
}

/*
 * (see inherited doc)
 */
ConfigDatabasePtr MySqlConfigDao::searchDatabase(const std::string& dbName)
{
	Statement stmt = (session->prepare <<
		"SELECT schema_name "
		"FROM information_schema.schemata "
		"WHERE schema_name=:name",
		use(const_cast<std::string&>(dbName)));
	return loadDatabase(stmt);
}

/*
 * (see inherited doc)
 */
ConfigDatabasePtr MySqlConfigDao::getCurrentDatabase()
{
	ConfigDatabasePtr dbObj(new ConfigDatabase);
	dbObj->setDAO(this);
	session->once << "SELECT database() AS currentDb", into(dbObj->name);
	return dbObj;
}

/*
 * (see inherited doc)
 */
ConfigDatabasePtrSetPtr MySqlConfigDao::getDatabases()
{
	ConfigDatabasePtrSetPtr databases(new ConfigDatabasePtrSet);
	Rowset<ConfigDatabasePtr> rowset =
		(session->prepare << "SELECT schema_name FROM information_schema.schemata");
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
std::vector<std::string> MySqlConfigDao::getUserNames()
{
	std::vector<std::string> users;
	Rowset<std::string> rowset = (session->prepare <<
			"SELECT User FROM mysql.user WHERE Host='%'");
	for(Rowset<std::string>::iterator it = rowset.begin(); it != rowset.end(); ++it)
	{
		users.push_back(*it);
	}
	return users;
}

/*
 * (see inherited doc)
 */
void MySqlConfigDao::createDatabaseAndUser(const std::string& name, const std::string& pwd)
{
	session->once << "CREATE DATABASE " << name;
	session->once << "CREATE USER " << name << " IDENTIFIED BY :pwd ", use(const_cast<std::string&>(pwd));
	session->once << "GRANT ALL ON " << name << ".* TO " << name << " WITH GRANT OPTION";
	session->once << "GRANT SELECT(Host,User) ON mysql.user TO " << name;
	session->once << "GRANT SELECT ON mysql.db TO " << name;
	session->once << "GRANT SELECT ON mysql.tables_priv TO " << name;
}

/*
 * (see inherited doc)
 */
void MySqlConfigDao::removeDatabaseAndUser(const std::string& name)
{
	session->once << "REVOKE ALL PRIVILEGES, GRANT OPTION FROM " << name; // in MySQL only users without any rights can be dropped
	session->once << "DROP USER " << name;
	session->once << "DROP DATABASE " << name;
}

/*
 * (see inherited doc)
 */
void MySqlConfigDao::changeUserPassword(const std::string& userName, const std::string& newPassword)
{
	session->once <<
		"SET PASSWORD FOR " << userName << " = PASSWORD(:pwd)",
		use(const_cast<std::string&>(newPassword));
}

/*
 * (see inherited doc)
 */
int MySqlConfigDao::getPrivileges(const std::string& userName, const std::string& tableName)
{
	eIndicator ind;
	int privileges = 0;
	session->once <<
		"SELECT (LOCATE('Select',Table_priv)<>0) | "
			"((LOCATE('Insert',Table_priv)<>0)<<1) | "
			"((LOCATE('Update',Table_priv)<>0)<<2) | "
			"((LOCATE('Delete',Table_priv)<>0)<<3) AS Flags "
		"FROM mysql.tables_priv "
		"WHERE table_name=:table AND Db=:schema AND User=:user AND host='%'",
		use(const_cast<std::string&>(tableName)),
		use(const_cast<std::string&>(dbName)),
		use(const_cast<std::string&>(userName)),
		into(privileges, ind);
	return privileges;
}

/*
 * (see inherited doc)
 */
void MySqlConfigDao::revokePrivileges(const std::string& userName)
{
	for(size_t i = 0; i < TableCount+1; ++i)
	{
		session->once <<
			"REVOKE ALL PRIVILEGES, GRANT OPTION "
			"ON " << dbName << '.' << TABLES[i] << " "
			"FROM " << userName;
	}
}

/*
 * (see inherited doc)
 */
bool MySqlConfigDao::isSchemaAdmin(const std::string& userName)
{
	eIndicator ind;
	int isAdmin = 0;
	session->once <<
		"SELECT Select_priv = 'Y' AND "
			"Insert_priv = 'Y' AND "
			"Update_priv = 'Y' AND "
			"Delete_priv = 'Y' "
		"FROM mysql.db "
		"WHERE Db=:schema AND User=:user AND host='%'",
		use(const_cast<std::string&>(dbName)),
		use(const_cast<std::string&>(userName)),
		into(isAdmin, ind);
	return isAdmin != 0;
}

/*
 * (see inherited doc)
 */
bool MySqlConfigDao::checkTableExist(const char* tableName)
{
	int exist = -1;
	std::string name(tableName);
	session->once << "SELECT COUNT(*) AS tableExists "
		"FROM information_schema.TABLES "
		"WHERE TABLE_SCHEMA=:schema AND TABLE_NAME=:name",
		use(dbName),
		use(name),
		into(exist);
	return exist == 1;
}

/*
 * (see inherited doc)
 */
std::auto_ptr<TableBuilder> MySqlConfigDao::createTableBuilder()
{
	return std::auto_ptr<TableBuilder>(new MySqlTableBuilder);
}

}
}
}
