#include "Sqlite3ConfigDao.h"
#include "Sqlite3TableBuilder.h"
#include "ObjRelMapping.hpp"
#include "Transaction.h"
#include <sqlite3/soci-sqlite3.h>
#include <soci.h>
#include <sstream>
#include <cassert>

namespace CCfg
{
namespace Io
{
namespace Db
{
using namespace SOCI;

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::Sqlite3ConfigDao class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
Sqlite3ConfigDao::Sqlite3ConfigDao() : ConfigDao()
{
}

/**
 * Destructor.
 */
Sqlite3ConfigDao::~Sqlite3ConfigDao()
{
}

/*
 * (see inherited doc)
 */
void Sqlite3ConfigDao::connect(const std::string& dbPath)
{
	session.reset(new Session(sqlite3, dbPath));
	dbName = "main";
}

/*
 * (see inherited doc)
 */
ConfigDatabasePtr Sqlite3ConfigDao::searchDatabase(const std::string& dbName)
{
	ConfigDatabasePtr database(0);
	if (dbName == "main")
		database = getCurrentDatabase();
	return database;
}

/*
 * (see inherited doc)
 */
ConfigDatabasePtr Sqlite3ConfigDao::getCurrentDatabase()
{
	ConfigDatabasePtr dbObj(new ConfigDatabase);
	dbObj->setDAO(this);
	dbObj->name = "main";
	return dbObj;
}

/*
 * (see inherited doc)
 */
ConfigDatabasePtrSetPtr Sqlite3ConfigDao::getDatabases()
{
	ConfigDatabasePtrSetPtr databases(new ConfigDatabasePtrSet);
	databases->insert(getCurrentDatabase());
	return databases;
}

/*
 * (see inherited doc)
 */
std::vector<std::string> Sqlite3ConfigDao::getUserNames()
{
	assert("User management is not supported by SQLite 3.");
	return std::vector<std::string>();
}

/*
 * (see inherited doc)
 */
void Sqlite3ConfigDao::createDatabaseAndUser(const std::string&, const std::string&)
{
	assert("User management is not supported by SQLite 3.");
}

/*
 * (see inherited doc)
 */
void Sqlite3ConfigDao::removeDatabaseAndUser(const std::string&)
{
	assert("User management is not supported by SQLite 3.");
}

/*
 * (see inherited doc)
 */
void Sqlite3ConfigDao::changeUserPassword(const std::string&, const std::string&)
{
	assert("User management is not supported by SQLite 3.");
}

/*
 * (see inherited doc)
 */
int Sqlite3ConfigDao::getPrivileges(const std::string&, const std::string&)
{
	const int privileges = Select | Insert | Update | Delete;
	return privileges;
}

/*
 * (see inherited doc)
 */
void Sqlite3ConfigDao::revokePrivileges(const std::string&)
{
	;
}

/*
 * (see inherited doc)
 */
bool Sqlite3ConfigDao::isSchemaAdmin(const std::string&)
{
	return true;
}

/*
 * (see inherited doc)
 */
bool Sqlite3ConfigDao::checkTableExist(const char* tableName)
{
	int exist = -1;
	std::string name(tableName);
	session->once << "SELECT COUNT(*) AS tableExists "
		"FROM " << dbName << ".sqlite_master "
		"WHERE name=:name",
		use(name),
		into(exist);
	return exist == 1;
}

/**
 * Re-implements the allocation of Unique Identifiers for SQLite3, which does not support "SELECT FOR UPDATE" statements.
 * See inherited documentation.
 */
UID Sqlite3ConfigDao::allocIdRange(const char* tableName, size_t count)
{
	UID rangeBegin;
	Transaction transac(*session); // Should we use IMMEDIATE transactions?
	session->once <<
		"SELECT " << tableName << " "
		"FROM " << dbName << '.' << SEQUENCE_TABLE_NAME,
		into(rangeBegin);
	session->once <<
		"UPDATE " << dbName << '.' << SEQUENCE_TABLE_NAME << " "
		"SET " << tableName << '=' << tableName << '+' << count;
	transac.commit();
	return rangeBegin;
}

/**
 * Re-implements the allocation of Unique Identifiers for SQLite3, which does not support "SELECT FOR UPDATE" statements.
 * See inherited documentation.
 */
void Sqlite3ConfigDao::allocIdRanges(const UID (&counts)[TableCount], UID (&ids)[TableCount])
{
	Transaction transac(*session); // Should we use IMMEDIATE transactions?
	session->once <<
		"SELECT " <<
			TABLES[TConfig] << ',' <<
			TABLES[TStructureVersion] << ',' <<
			TABLES[TVersion] << ',' <<
			TABLES[TClassInfo] << ',' <<
			TABLES[TClass] << ',' <<
			TABLES[TAttribute] << ',' <<
			TABLES[TSite] << ',' <<
			TABLES[TRun] << ',' <<
			TABLES[TRunCondition] << " "
		"FROM " << dbName << '.' << SEQUENCE_TABLE_NAME,
		into(ids[TConfig]),
		into(ids[TStructureVersion]),
		into(ids[TVersion]),
		into(ids[TClassInfo]),
		into(ids[TClass]),
		into(ids[TAttribute]),
		into(ids[TSite]),
		into(ids[TRun]),
		into(ids[TRunCondition]);
	session->once <<
		"UPDATE " << dbName << '.' << SEQUENCE_TABLE_NAME << " "
		"SET " <<
			TABLES[TConfig] << '=' << TABLES[TConfig] << '+' << counts[TConfig] << ',' <<
			TABLES[TStructureVersion] << '=' << TABLES[TStructureVersion] << '+' << counts[TStructureVersion] << ',' <<
			TABLES[TVersion] << '=' << TABLES[TVersion] << '+' << counts[TVersion] << ',' <<
			TABLES[TClassInfo] << '=' << TABLES[TClassInfo] << '+' << counts[TClassInfo] << ',' <<
			TABLES[TClass] << '=' << TABLES[TClass] << '+' << counts[TClass] << ',' <<
			TABLES[TAttribute] << '=' << TABLES[TAttribute] << '+' << counts[TAttribute] << ',' <<
			TABLES[TSite] << '=' << TABLES[TSite] << '+' << counts[TSite] << ',' <<
			TABLES[TRun] << '=' << TABLES[TRun] << '+' << counts[TRun] << ',' <<
			TABLES[TRunCondition] << '=' << TABLES[TRunCondition] << '+' << counts[TRunCondition];
	transac.commit();
}

/*
 * (see inherited doc)
 */
std::auto_ptr<TableBuilder> Sqlite3ConfigDao::createTableBuilder()
{
	return std::auto_ptr<TableBuilder>(new Sqlite3TableBuilder);
}

}
}
}
